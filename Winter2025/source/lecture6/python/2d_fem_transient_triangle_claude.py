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
from fem_utils import compute_jacobian
from mesh_utils import make_box, make_inclusion
from validation import check_cfl_condition, validate_parameters

import triangle as tr
import meshio 

# Model parameters
#geometry
x0          = -1
y0          = -1
lx          = 2
ly          = 2
n_incl      = 5
radius      = 0.15

# time control
dt = 0.015  #arbitrary time stepping
nt = 80

# model parameters
rho         = 1
cp          = 1
k1          = 1
k2          = 0.01
Ttop        = 0
Tbot        = 1

# Validate parameters
params = {
    'rho': rho,
    'cp': cp,
    'k1': k1,
    'k2': k2,
    'dt': dt
}
validate_parameters(params)

## Create the triangle mesh
# arrays to fill in with input
vertices = []
segments = []
regions = []

# make_box and make_inclusion are now imported from mesh_utils (removed inline definitions)

# generate input - using utility functions
make_box(vertices, segments, regions, x0, y0, lx, ly, attribute=1, max_area=0.005)

make_inclusion(vertices, segments, regions, -0.8, -0.3, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, -0.5, -0.75, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, -0.6, 0.5, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, -0.1, -0.3, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, 0.1, 0, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, 0.5, -0.2, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, 0.6, .3, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, 0.7, .8, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, 0, .75, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, -0.5, .05, radius, 20, attribute=100, max_area=0.001)
make_inclusion(vertices, segments, regions, 0.5, -.75, radius, 20, attribute=100, max_area=0.001)

A = dict(vertices=vertices, segments=segments, regions=regions)
B = tr.triangulate(A, 'pq33Aa')

# extract mesh information
GCOORD = B.get("vertices")
EL2NOD = B.get("triangles")
Phases = B.get("triangle_attributes")

nnodel = EL2NOD.shape[1]
nel    = EL2NOD.shape[0]
nnod   = GCOORD.shape[0]
Phases = np.reshape(Phases,nel)
Kel    = np.ones(nel)*k1
Kel[np.where(Phases==100)] = k2
T           = np.zeros(nnod) #initial T

# Estimate characteristic mesh spacing for stability check
coords_x = GCOORD[:, 0]
coords_y = GCOORD[:, 1]
dx = np.sqrt((np.max(coords_x) - np.min(coords_x)) * (np.max(coords_y) - np.min(coords_y)) / nnod)

# Check CFL condition (informational - implicit method is unconditionally stable)
# Use the larger thermal diffusivity for the check
alpha_max = max(k1/(rho*cp), k2/(rho*cp))
is_stable, cfl, cfl_max = check_cfl_condition(dt, dx, alpha_max, method='implicit', ndim=2)
print(f"CFL number: {cfl:.4f} (max for explicit: {cfl_max:.4f})")
print(f"Using implicit method - unconditionally stable")

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

    # Storage
    Rhs_all = np.zeros(nnod)
    I       = np.zeros((nel,nnodel*nnodel))
    J       = np.zeros((nel,nnodel*nnodel))
    K       = np.zeros((nel,nnodel*nnodel))

    for iel in range(0,nel):
        ECOORD = np.take(GCOORD, EL2NOD[iel,:], axis=0 )
        Ael    = np.zeros((nnodel,nnodel))
        Rhs_el = np.zeros(nnodel)
        
        for ip in range(0,nip):
            # 1. update shape functions
            xi      = gauss[ip,0]
            eta     = gauss[ip,1]
            N, dNds = shapes_tri(xi, eta, nnodel)

            # 2. set up Jacobian, inverse of Jacobian, and determinant
            Jac, invJ, detJ = compute_jacobian(dNds, ECOORD)
            
            # 3. get global derivatives
            dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]
            
            # 4. compute element stiffness matrix
            Ael     = Ael + (rho*cp*np.outer(N,N) +  dt*Kel[iel]*np.matmul(dNdx.T, dNdx))*detJ*weights[ip] # [nnodel,1]*[1,nnodel] / weights are missing, they are 1
            
            # 5. assemble right-hand side
            Rhs_el     = Rhs_el + rho*cp*np.matmul(np.outer(N,N), np.take(T, EL2NOD[iel,:], axis=0 ))*detJ*weights[ip] 
        
        # assemble coefficients
        I[iel,:]  =  (EL2NOD[iel,:]*np.ones((nnodel,1), dtype=int)).T.reshape(nnodel*nnodel)
        J[iel,:]  =  (EL2NOD[iel,:]*np.ones((nnodel,1), dtype=int)).reshape(nnodel*nnodel)
        K[iel,:]  =  Ael.reshape(nnodel*nnodel)
        
        Rhs_all[EL2NOD[iel,:]] += Rhs_el

    A_all = csr_matrix((K.reshape(nel*nnodel*nnodel),(I.reshape(nel*nnodel*nnodel),J.reshape(nel*nnodel*nnodel))),shape=(nnod,nnod))

    # indices and values at top and bottom
    tol = 1e-3
    i_bot = np.where(abs(GCOORD[:,1] - y0) < tol)[0]
    i_top = np.where(abs(GCOORD[:,1] - (y0+lx)) < tol)[0]

    Ind_bc  = np.concatenate((i_bot, i_top))
    Val_bc  = np.concatenate((np.ones(i_bot.shape)*Tbot, np.ones(i_top.shape)*Ttop ))

    # boundary conditions that keep matrix symmetry
    Free    = np.arange(0,nnod)
    Free    = np.delete(Free, Ind_bc)
    TMP     = A_all[:,Ind_bc]
    Rhs_all = Rhs_all - TMP.dot(Val_bc)

    # solve reduced system
    T[Free] = spsolve(A_all[np.ix_(Free, Free)],Rhs_all[Free])
    T[Ind_bc] = Val_bc

    # postprocessing - heat flow
    Q_x  = np.zeros(nel)
    Q_y  = np.zeros(nel)
    Ec_x = np.zeros(nel)
    Ec_y = np.zeros(nel)

    for iel in range(0,nel):
        # 0. get element coordinates
        ECOORD = np.take(GCOORD, EL2NOD[iel,:], axis=0 )        
        # 1. update shape functions
        xi      = 1/3
        eta     = 1/3
        N, dNds = shapes_tri(xi, eta)           
        # 2. set up Jacobian, inverse of Jacobian, and determinant
        Jac     = np.matmul(dNds,ECOORD) #[2,nnodel]*[nnodel,2]
        invJ    = np.linalg.inv(Jac)    
        detJ    = np.linalg.det(Jac)           
        # 3. get global derivatives
        dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]           
        # 4. heat flux per element  
        Q_x[iel] = -Kel[iel]*np.matmul(dNdx[0,:], np.take(T, EL2NOD[iel,:]))
        Q_y[iel] = -Kel[iel]*np.matmul(dNdx[1,:], np.take(T, EL2NOD[iel,:]))
        Ec_x[iel]  = np.mean(ECOORD[:,0])
        Ec_y[iel]  = np.mean(ECOORD[:,1])

    #cell data
    U=np.hstack((Q_x.reshape(-1,1),Q_y.reshape(-1,1)))
    U=np.hstack((U,U[:,0].reshape(-1,1)*0))
    
    #save data
    writer.write_data(t, point_data={"T": T},cell_data={"U": [U], "K": [Kel]})

writer.__exit__() # close file



