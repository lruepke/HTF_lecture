import numpy as np
from shp_deriv_triangle import *
from ip_triangle import *
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def mechanical2d(Mesh, Materials, Solver):
    """
    MECHANICAL2D Two dimensional finite element mechanical problem solver.

    Args:
        Mesh (dict): Mesh parameters containing:
            - 'EL2NOD': Element-to-node connectivity.
            - 'Phases': Phases.
            - 'GCOORD': Global coordinates.
            - 'Bc_ind': Boundary condition indices.
            - 'Bc_val': Boundary condition values.
        Materials (dict): Material properties containing:
            - 'Mu': Shear modulus.
            - 'Rho': Density.
            - 'G': Gravity vector.
        Solver (dict): Solver switches:
            - 'nip': number of integration points.
    """
    
    # Unpack mesh parameters
    EL2NOD  = Mesh.get('EL2NOD')
    Phases  = Mesh.get('Phases')
    GCOORD  = Mesh.get('GCOORD')
    Bc_ind  = Mesh.get('Bc_ind')
    Bc_val  = Mesh.get('Bc_val')
    nel     = Mesh.get('nel')
    nnod    = Mesh.get('nnod')
    nnodel  = Mesh.get('nnodel')
    
    # Unpack material properties
    Mu      = Materials['D']
    Rho     = Materials['Rho']
    G       = Materials['G']
    
    # Unpack solver settings
    nip     = Solver['nip']
    
    # CONSTANTS
    ndim = 2
    nedof = nnodel * ndim
    sdof = 2 * nnod
    np_edof = 3 #number of pressurse eldof

    DEV = np.array([[4/3, -2/3, 0],
                    [-2/3, 4/3, 0],
                    [0, 0, 1]])

    PF = 1e3 * np.max(Mu)

    C1 = 4/3
    C2 = 2/3

    # PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
    IP_X, IP_w = ip_triangle(nip)
    Shapes = shp_deriv_triangle(IP_X, nnodel)

    # DECLARE VARIABLES (ALLOCATE MEMORY)
    A_all   = np.zeros((nel, nedof*nedof))
    I       = np.zeros((nel, nedof*nedof))
    J       = np.zeros((nel, nedof*nedof))
    Q_all   = np.zeros((nel, nedof * np_edof))
    invM_all= np.zeros((nel,np_edof*np_edof)) 
    Rhs_all = np.zeros(sdof)
    Vel     = np.zeros(sdof) #initial T, not strictly needed

    # MAKE ELDOF
    EL2DOF = np.zeros((nel,nedof), dtype=np.int32)
    EL2DOF[:,0::ndim] = ndim * EL2NOD
    EL2DOF[:,1::ndim] = ndim * EL2NOD + 1

    # ELEMENT LOOP - MATRIX COMPUTATION
    for iel in range(nel):

        # FETCH DATA OF ELEMENT
        ECOORD_X = GCOORD[EL2NOD[iel,:],:]
        EMu = Mu[Phases[iel]]
        ERho = Rho[Phases[iel]]

        # INTEGRATION LOOP
        A_elem = np.zeros((nedof, nedof))
        Q_elem = np.zeros((nedof, np_edof))
        M_elem = np.zeros((np_edof, np_edof))
        Rhs_elem = np.zeros(nedof)

        B = np.zeros((ndim*(ndim+1)//2, nedof))
        P = np.ones((np_edof, np_edof))
        Pb = np.ones(np_edof)

        P[:, 1:3] = ECOORD_X[:3]
        for ip in range(nip):
            # LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
            Ni      = Shapes['shape_functions'][ip]
            dNdui   = Shapes['shape_func_deriv'][ip]
            Pb[1:3] = Ni @ ECOORD_X
            Pi      = np.linalg.solve(P, Pb)

            # CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
            Jac       = dNdui @ ECOORD_X

            detJ    = np.linalg.det(Jac)
            invJ    = np.linalg.inv(Jac)

            # DERIVATIVES wrt GLOBAL COORDINATES
            dNdX    = invJ @ dNdui 

            # NUMERICAL INTEGRATION OF ELEMENT MATRICES
            weight = IP_w[ip] * detJ
            B[0,0::2] = dNdX[0]
            B[1,1::2] = dNdX[1]
            B[2,0::2] = dNdX[1]
            B[2,1::2] = dNdX[0]
            m         = np.array([1,1,0])

            A_elem += weight * EMu * (B.T @ DEV @ B)
            Q_elem -= weight * B.T @ m[:, np.newaxis] @ Pi[np.newaxis, :]
            M_elem += weight * Pi[:, np.newaxis] @ Pi[np.newaxis, :]
            Rhs_elem += weight * ERho * (G[:, np.newaxis] @ Ni[np.newaxis, :]).T.ravel()

        # STATIC CONDENSATION
        invM_elem = np.linalg.inv(M_elem)
        A_elem += PF * (Q_elem @ invM_elem @ Q_elem.T)

        # WRITE DATA INTO GLOBAL STORAGE
        A_all[iel, :] = A_elem.ravel() #A_elem[np.triu_indices(nedof)]
        Q_all[iel,:] = Q_elem.ravel()
        invM_all[iel,:] = invM_elem.ravel()
        Rhs_all[EL2DOF[iel,:]] += Rhs_elem

        # assemble coefficients
        I[iel,:]  =  (EL2DOF[iel,:]*np.ones((nedof,1), dtype=int)).T.reshape(nedof*nedof)
        J[iel,:]  =  (EL2DOF[iel,:]*np.ones((nedof,1), dtype=int)).reshape(nedof*nedof)


    # ASSEMBLE GLOBAL MATRICES
    A_all = csr_matrix((A_all.reshape(nel*nedof*nedof),(I.reshape(nel*nedof*nedof),J.reshape(nel*nedof*nedof))),shape=(sdof,sdof))


    # smart way of boundary conditions that keeps matrix symmetry
    Free    = np.arange(0,sdof)
    Free    = np.delete(Free, Bc_ind)
    TMP     = A_all[:,Bc_ind]
    Rhs_all = Rhs_all - TMP.dot(Bc_val)

    # solve reduced system
    Vel[Free] = spsolve(A_all[np.ix_(Free, Free)],Rhs_all[Free])
    Vel[Bc_ind] = Bc_val

    # compute pressure
    

    return Vel