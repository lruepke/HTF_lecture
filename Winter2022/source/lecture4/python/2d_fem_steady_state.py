import numpy as np
from tabulate import tabulate
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
import matplotlib as mpl
import matplotlib.pyplot as plt
from shapes import shapes

#geometry
lx          = 1
ly          = 1
nx          = 51
ny          = 51
nnodel      = 4
dx          = lx/(nx-1)
dy          = ly/(ny-1)
w_i         = 0.2 # width inclusion
h_i         = 0.2 # heigths inclusion

# model parameters
k1          = 1
k2          = 0.001
Ttop        = 0
Tbot        = 1
 
nex         = nx-1
ney         = ny-1
nnod        = nx*ny
nel         = nex*ney
GCOORD      = np.zeros((nnod,2))
T           = np.zeros(nnod) #initial T, not strictly needed

# global coordinates
id = 0
for i in range(0,ny):
    for j in range(0,nx):
        GCOORD[id,0] = -lx/2 + j*dx
        GCOORD[id,1] = -ly/2 + i*dy
        id          = id + 1

# FEM connectivity
EL2NOD   = np.zeros((nel,nnodel), dtype=int)

for iel in range(0,nel):
    row        = iel//nex   
    ind        = iel + row
    EL2NOD[iel,:] = [ind, ind+1, ind+nx+1, ind+nx]
    
# Gauss integration points
nip   = 4
gauss = np.array([[ -np.sqrt(1/3), np.sqrt(1/3), np.sqrt(1/3), -np.sqrt(1/3)], [-np.sqrt(1/3), -np.sqrt(1/3), np.sqrt(1/3), np.sqrt(1/3)]]).T.copy()

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
        N, dNds = shapes(xi, eta)
        
        # 2. set up Jacobian, inverse of Jacobian, and determinant
        Jac     = np.matmul(dNds,ECOORD) #[2,nnodel]*[nnodel,2]
        invJ    = np.linalg.inv(Jac)     
        detJ    = np.linalg.det(Jac)
        
        # 3. get global derivatives
        dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]
        
        # 4. compute element stiffness matrix
        kel = k1
        if  abs(np.mean(ECOORD[:,0]))<w_i and abs(np.mean(ECOORD[:,1]))<h_i:
            kel=k2  
        Ael     = Ael + np.matmul(dNdx.T, dNdx)*detJ*kel # [nnodel,1]*[1,nnodel] / weights are missing, they are 1
        
        # 5. assemble right-hand side, no source terms, just here for completeness
        Rhs_el     = Rhs_el + np.zeros(nnodel)
    
    # assemble coefficients
    I[iel,:]  =  (EL2NOD[iel,:]*np.ones((nnodel,1), dtype=int)).T.reshape(nnodel*nnodel)
    J[iel,:]  =  (EL2NOD[iel,:]*np.ones((nnodel,1), dtype=int)).reshape(nnodel*nnodel)
    K[iel,:]  =  Ael.reshape(nnodel*nnodel)
    
    Rhs_all[EL2NOD[iel,:]] += Rhs_el

A_all = csr_matrix((K.reshape(nel*nnodel*nnodel),(I.reshape(nel*nnodel*nnodel),J.reshape(nel*nnodel*nnodel))),shape=(nnod,nnod))

# indices and values at top and bottom
i_bot   = np.arange(0,nx, dtype=int)
i_top   = np.arange(nx*(ny-1),nx*ny, dtype=int)
Ind_bc  = np.concatenate((i_bot, i_top))
Val_bc  = np.concatenate((np.ones(i_bot.shape)*Tbot, np.ones(i_top.shape)*Ttop ))

# smart way of boundary conditions that keeps matrix symmetry
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
    xi      = 0
    eta     = 0
    N, dNds = shapes(xi, eta)
    
    # 2. set up Jacobian, inverse of Jacobian, and determinant
    Jac     = np.matmul(dNds,ECOORD) #[2,nnodel]*[nnodel,2]
    invJ    = np.linalg.inv(Jac)    
    detJ    = np.linalg.det(Jac)
    
    # 3. get global derivatives
    dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]
    
    # 4. heat flux per element
    kel = k1
    if  abs(np.mean(ECOORD[:,0]))<0.25 and abs(np.mean(ECOORD[:,1]))<0.25:
        kel=k2  
    Q_x[iel] = -kel*np.matmul(dNdx[0,:], np.take(T, EL2NOD[iel,:]))
    Q_y[iel] = -kel*np.matmul(dNdx[1,:], np.take(T, EL2NOD[iel,:]))
    Ec_x[iel]  = np.mean(ECOORD[:,0])
    Ec_y[iel]  = np.mean(ECOORD[:,1])

# plotting
fig = plt.figure()
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height]) 

X = np.reshape(GCOORD[:,0], (ny,nx))
Y = np.reshape(GCOORD[:,1], (ny,nx))
T = T.reshape((ny,nx))

cp = plt.contourf(X, Y, T, cmap='gist_yarg')
plt.colorbar(cp)
plt.quiver(Ec_x, Ec_y, Q_x, Q_y, np.sqrt(np.square(Q_x) + np.square(Q_y)), cmap='hot')

ax.set_title('Temperature with heat flow vectors')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()



