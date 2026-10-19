# Script that calculated steady-state diffusion on unstructured triangle mesh
import numpy as np
from tabulate import tabulate
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
import matplotlib as mpl
import matplotlib.pyplot as plt

# FEM utilities (refactored 2025)
from fem_shapes import shapes_tri
from fem_integration import ip_triangle
from fem_utils import compute_jacobian, create_assembly_indices_2dof, assemble_sparse_matrix_2dof
from mesh_utils import make_box
from validation import check_cfl_condition, validate_parameters

import triangle as tr
import meshio
from numpy.random import default_rng

# Model parameters
#geometry
x0          = 0
y0          = 0
lx          = 5
ly          = 5

# time control
dt = 5e-4
nt = 200

# model parameters
g_coeff     = 600 #kinetic constant
a_coeff     = 0.05 #growth A
b_coeff     = 1    #growth B
d_coeff     = 20   #diffusivity B
amp         = 0.01 #random noise

# Validate parameters
params = {
    'd_coeff': d_coeff,
    'dt': dt,
    'g_coeff': g_coeff,
    'a_coeff': a_coeff,
    'b_coeff': b_coeff
}
validate_parameters(params)

## Create the triangle mesh
# arrays to fill in with input
vertices = []
segments = []
regions = []

# make_box is now imported from mesh_utils (removed inline definition)

# generate input - using utility function
make_box(vertices, segments, regions, x0, y0, lx, ly, attribute=1, max_area=0.005)
A = dict(vertices=vertices, segments=segments, regions=regions)
B = tr.triangulate(A, 'pq33Aa')

# extract mesh information
GCOORD = B.get("vertices")
EL2NOD = B.get("triangles")
Phases = B.get("triangle_attributes")

nnodel = EL2NOD.shape[1]
nel    = EL2NOD.shape[0]
nnod   = GCOORD.shape[0]
sdof   = nnod*2                 # two dof per node
Phases = np.reshape(Phases,nel)
print(nnod, nel)

# Estimate characteristic mesh spacing for stability check
coords_x = GCOORD[:, 0]
coords_y = GCOORD[:, 1]
dx = np.sqrt((np.max(coords_x) - np.min(coords_x)) * (np.max(coords_y) - np.min(coords_y)) / nnod)

# Check CFL condition (informational - implicit method is unconditionally stable)
is_stable, cfl, cfl_max = check_cfl_condition(dt, dx, d_coeff, method='implicit', ndim=2)
print(f"CFL number: {cfl:.4f} (max for explicit: {cfl_max:.4f})")
print(f"Using implicit method - unconditionally stable")
# setup degrees of freedom - two per node
EL2DOF = np.zeros((nel,2*nnodel), dtype=int)
EL2DOF[:,0::2] = 2*EL2NOD
EL2DOF[:,1::2] = 2*EL2NOD+1

# Initial conditions
rng  = default_rng()
vals = rng.standard_normal(nnod)
A    = (a_coeff+b_coeff) + amp*vals
vals = rng.standard_normal(nnod)
B    = b_coeff/(a_coeff+b_coeff)**2 + amp*vals

# setup output writing
points=np.hstack((GCOORD, GCOORD[:,0].reshape(-1,1)*0)) #must have 3 components (x,y,z)
cells=[("triangle",EL2NOD)]
writer=meshio.xdmf.TimeSeriesWriter('transient.xmf')
writer.__enter__() # have to add this: import hdf5 and open file ...
writer.write_points_cells(points, cells)

# Gauss integration points for triangles
nip   = 3
gauss, weights = ip_triangle(nip)

# time loop

for t in range(0,nt):

    # Storage - use utility function for assembly indices
    Rhs_all = np.zeros(sdof)
    I, J = create_assembly_indices_2dof(EL2NOD, nnodel)
    K    = np.zeros((nel,2*nnodel*nnodel))

    for iel in range(0,nel):
        ECOORD  = np.take(GCOORD, EL2NOD[iel,:], axis=0 )
        Ael_A   = np.zeros((nnodel,nnodel))
        Ael_B   = np.zeros((nnodel,nnodel))
        RhsA_el = np.zeros(nnodel)
        RhsB_el = np.zeros(nnodel)
        
        for ip in range(0,nip):        
            # 1. update shape functions
            xi      = gauss[ip,0]
            eta     = gauss[ip,1]
            N, dNds = shapes_tri(xi, eta, nnodel)

            # 2. set up Jacobian, inverse of Jacobian, and determinant
            Jac, invJ, detJ = compute_jacobian(dNds, ECOORD)
            
            # 3. get global derivatives
            dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]
            
            # 4. compute element stiffness matrices
            Ael_A     = Ael_A + (np.outer(N,N)*(1+g_coeff*dt) +  dt*np.matmul(dNdx.T, dNdx))*detJ*weights[ip] 
            Ael_B     = Ael_B + (np.outer(N,N)        +  d_coeff*dt*np.matmul(dNdx.T, dNdx))*detJ*weights[ip] 
            
            # 5. assemble right-hand side
            RhsA_el     = RhsA_el + np.matmul(np.outer(N,N), np.take(A, EL2NOD[iel,:], axis=0 ))*detJ*weights[ip] 
            RhsB_el     = RhsB_el + np.matmul(np.outer(N,N), np.take(B, EL2NOD[iel,:], axis=0 ))*detJ*weights[ip] 
        
        # assemble element stiffness coefficients (I, J already created by utility function)
        K[iel,:]  =  np.concatenate((Ael_A.reshape(nnodel*nnodel),Ael_B.reshape(nnodel*nnodel)))
        
        Rhs_all[2*EL2NOD[iel,:]]   += RhsA_el
        Rhs_all[2*EL2NOD[iel,:]+1] += RhsB_el

    A_all = csr_matrix((K.reshape(nel*2*nnodel*nnodel),(I.reshape(nel*2*nnodel*nnodel),J.reshape(nel*2*nnodel*nnodel))),shape=(sdof,sdof))

    
    # update right hand side in iterations

    iter = 0
    error = 10
    tol   = 0.001
    Conc_tmp = np.ones(sdof)*10
    iter_max = 20

    while error > tol:
        Tmp = Rhs_all.copy()
        iter += 1
        # loop over all elements and integrate Rhs
        for iel in range(0,nel):
            FA_el = np.zeros(nnodel)
            FB_el = np.zeros(nnodel)
            ECOORD  = np.take(GCOORD, EL2NOD[iel,:], axis=0 )

            for ip in range(0,nip):
                # 1. update shape functions
                xi      = gauss[ip,0]
                eta     = gauss[ip,1]
                N, dNds = shapes_tri(xi, eta, nnodel)

                # 2. set up Jacobian, inverse of Jacobian, and determinant
                Jac, invJ, detJ = compute_jacobian(dNds, ECOORD)

                #3. integrate force vector
                Ai = np.dot(N,np.take(A, EL2NOD[iel,:], axis=0 ))
                Bi = np.dot(N,np.take(B, EL2NOD[iel,:], axis=0 ))
                FA_el     = FA_el + N*dt*g_coeff*(a_coeff+Bi*Ai**2)*detJ*weights[ip] # (dt*g_coeff*N*a_coeff+dt*g_coeff*N*np.dot(N,np.take(A, EL2NOD[iel,:], axis=0 ))**2*np.dot(N,np.take(B, EL2NOD[iel,:], axis=0 )))*detJ*weights[ip] 
                FB_el     = FB_el + N*dt*g_coeff*(b_coeff-Bi*Ai**2)*detJ*weights[ip] # (dt*g_coeff*N*b_coeff-dt*g_coeff*N*np.dot(N,np.take(A, EL2NOD[iel,:], axis=0 ))**2*np.dot(N,np.take(B, EL2NOD[iel,:], axis=0 )))*detJ*weights[ip] 

                # 2. integrate force terms
      #          FA_el     = FA_el + (dt*g_coeff*N*a_coeff+dt*g_coeff*N*np.dot(N,np.take(A, EL2NOD[iel,:], axis=0 ))**2*np.dot(N,np.take(B, EL2NOD[iel,:], axis=0 )))*detJ*weights[ip] 
      #          FB_el     = FB_el + (dt*g_coeff*N*b_coeff-dt*g_coeff*N*np.dot(N,np.take(A, EL2NOD[iel,:], axis=0 ))**2*np.dot(N,np.take(B, EL2NOD[iel,:], axis=0 )))*detJ*weights[ip] 
      
            # We don't have boundary conditions, as everything is zero flux      
            Tmp[2*EL2NOD[iel,:]]   += FA_el
            Tmp[2*EL2NOD[iel,:]+1] += FB_el
       
        # solve  system
        Conc  = spsolve(A_all,Tmp)        
        error = np.amax(np.absolute(Conc - Conc_tmp))/np.amax(np.absolute(Conc))
        Conc_tmp = Conc.copy()
        A     = Conc[0:sdof:2]
        B     = Conc[1:sdof:2]
        
        print(error, iter)
        if iter == iter_max:

            break

      
    
    #save data
    print(t)
    writer.write_data(t, point_data={"A": A, "B": B})

writer.__exit__() # close file



