"""
2D Mechanical FEM solver for Stokes flow with discontinuous pressure.

This module solves 2D incompressible Stokes flow problems using triangular
finite elements with discontinuous pressure approximation.
"""

from dataclasses import dataclass
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from shp_deriv_triangle import shp_deriv_triangle
from ip_triangle import ip_triangle


@dataclass
class Mesh:
    """Container for mesh data."""
    GCOORD: np.ndarray   # Node coordinates, shape (nnod, 2)
    EL2NOD: np.ndarray   # Element connectivity, shape (nel, nnodel)
    Phases: np.ndarray   # Phase/material ID per element, shape (nel,)
    Node_ids: np.ndarray # Node boundary markers, shape (nnod,)
    nnod: int = 0        # Number of nodes
    nel: int = 0         # Number of elements
    nnodel: int = 0      # Number of nodes per element


@dataclass
class MaterialParams:
    """Container for material properties."""
    Mu: np.ndarray   # Shear viscosity per phase
    Rho: np.ndarray  # Density per phase
    G: np.ndarray    # Gravity vector [gx, gy]


@dataclass
class BoundaryConditions:
    """Container for boundary conditions."""
    Bc_ind: np.ndarray  # DOF indices with prescribed values
    Bc_val: np.ndarray  # Prescribed values at boundary DOFs


@dataclass
class SolverParams:
    """Container for solver parameters."""
    nip: int  # Number of integration points per element


@dataclass
class Solution:
    """Container for solution variables."""
    Vel: np.ndarray       # Velocity field, shape (2*nnod,)
    Pressure: np.ndarray  # Pressure field, shape (nel*3,)


def solve_mechanical_2d(
    mesh: Mesh,
    material: MaterialParams,
    bc: BoundaryConditions,
    solver: SolverParams
) -> Solution:
    """
    Solve 2D incompressible Stokes flow using FEM with discontinuous pressure.

    This solver uses triangular elements with quadratic velocity and linear
    discontinuous pressure (P2+/P-1 element). The incompressibility constraint
    is enforced using a penalty method with static condensation of pressure.

    Parameters
    ----------
    mesh : Mesh
        Mesh object containing coordinates, connectivity, and phase IDs
    material : MaterialParams
        Material properties (viscosity, density, gravity)
    bc : BoundaryConditions
        Boundary condition indices and values for velocity DOFs
    solver : SolverParams
        Solver parameters (number of integration points)

    Returns
    -------
    Solution
        Solution object containing velocity and pressure fields
    """
    # Extract mesh dimensions
    nel = mesh.nel
    nnod = mesh.nnod
    nnodel = mesh.nnodel
    
    # CONSTANTS
    ndim = 2
    nedof = nnodel * ndim      # DOFs per element
    sdof = 2 * nnod            # Total velocity DOFs
    np_edof = 3                # Pressure DOFs per element (linear on triangle)

    # Deviatoric projection operator for plane strain
    DEV = np.array([[4/3, -2/3, 0],
                    [-2/3, 4/3, 0],
                    [0, 0, 1]])

    # Penalty factor for incompressibility
    PF = 1e3 * np.max(material.Mu)

    # PREPARE INTEGRATION POINTS & SHAPE FUNCTIONS
    IP_X, IP_w = ip_triangle(solver.nip)
    Shapes = shp_deriv_triangle(IP_X, nnodel)

    # ALLOCATE STORAGE
    A_all = np.zeros((nel, nedof * nedof))
    I = np.zeros((nel, nedof * nedof))
    J = np.zeros((nel, nedof * nedof))
    Q_all = np.zeros((nel, nedof * np_edof))
    invM_all = np.zeros((nel, np_edof * np_edof))
    Rhs_all = np.zeros(sdof)

    # BUILD ELEMENT-TO-DOF MAPPING
    EL2DOF = np.zeros((nel, nedof), dtype=np.int32)
    EL2DOF[:, 0::ndim] = ndim * mesh.EL2NOD
    EL2DOF[:, 1::ndim] = ndim * mesh.EL2NOD + 1

    # ELEMENT ASSEMBLY LOOP
    for iel in range(nel):
        # Fetch element data
        ECOORD_X = mesh.GCOORD[mesh.EL2NOD[iel, :], :]
        EMu = material.Mu[mesh.Phases[iel]]
        ERho = material.Rho[mesh.Phases[iel]]

        # Initialize element matrices
        A_elem = np.zeros((nedof, nedof))
        Q_elem = np.zeros((np_edof, nedof))
        M_elem = np.zeros((np_edof, np_edof))
        Rhs_elem = np.zeros(nedof)

        # Strain-displacement matrix (B) and pressure interpolation
        B = np.zeros((ndim * (ndim + 1) // 2, nedof))
        P = np.ones((np_edof, np_edof))
        Pb = np.ones(np_edof)
        P[1:3, :] = ECOORD_X[:3].T

        # INTEGRATION POINT LOOP
        for ip in range(solver.nip):
            # Load shape functions and derivatives
            Ni = Shapes['shape_functions'][ip]
            dNdui = Shapes['shape_func_deriv'][ip]
            Pb[1:3] = Ni @ ECOORD_X
            Pi = np.linalg.solve(P, Pb)

            # Compute Jacobian and its inverse
            Jac = dNdui @ ECOORD_X
            detJ = np.linalg.det(Jac)
            invJ = np.linalg.inv(Jac)

            # Global derivatives of shape functions
            dNdX = invJ @ dNdui

            # Integration weight
            weight = IP_w[ip] * detJ

            # Build strain-displacement matrix B
            B[0, 0::2] = dNdX[0]  # dNdx for x-displacement
            B[1, 1::2] = dNdX[1]  # dNdy for y-displacement
            B[2, 0::2] = dNdX[1]  # dNdy for x-displacement (shear)
            B[2, 1::2] = dNdX[0]  # dNdx for y-displacement (shear)
            Bvol = dNdX.T         # Volumetric strain (divergence)

            # Accumulate element matrices
            A_elem += weight * EMu * (B.T @ DEV @ B)
            Q_elem -= weight * Pi[:, np.newaxis] @ Bvol.ravel()[np.newaxis, :]
            M_elem += weight * Pi[:, np.newaxis] @ Pi[np.newaxis, :]
            Rhs_elem += weight * ERho * (material.G[:, np.newaxis] @ Ni[np.newaxis, :]).T.ravel()

        # STATIC CONDENSATION of pressure DOFs
        invM_elem = np.linalg.inv(M_elem)
        A_elem += PF * (Q_elem.T @ invM_elem @ Q_elem)

        # Store element contributions
        A_all[iel, :] = A_elem.ravel()
        Q_all[iel, :] = Q_elem.ravel()
        invM_all[iel, :] = invM_elem.ravel()
        Rhs_all[EL2DOF[iel, :]] += Rhs_elem

        # Build sparse matrix indices
        I[iel, :] = np.tile(EL2DOF[iel, :], (nedof, 1)).T.ravel()
        J[iel, :] = np.tile(EL2DOF[iel, :], (nedof, 1)).ravel()

    # ASSEMBLE GLOBAL SPARSE MATRIX
    A_all = csr_matrix((A_all.ravel(), (I.ravel(), J.ravel())), shape=(sdof, sdof))

    # APPLY BOUNDARY CONDITIONS (preserving symmetry)
    Free = np.arange(0, sdof)
    Free = np.delete(Free, bc.Bc_ind)
    TMP = A_all[:, bc.Bc_ind]
    Rhs_all = Rhs_all - TMP.dot(bc.Bc_val)

    # SOLVE LINEAR SYSTEM
    Vel = np.zeros(sdof)
    Vel[Free] = spsolve(A_all[np.ix_(Free, Free)], Rhs_all[Free])
    Vel[bc.Bc_ind] = bc.Bc_val

    # RECOVER PRESSURE FIELD
    # Build global Q and invM matrices for discontinuous pressure
    Q_i = np.tile(np.arange(0, nel * np_edof, dtype=np.int32), (nedof, 1)).T
    Q_j = np.tile(EL2DOF, (1, np_edof))
    Q_all = csr_matrix((Q_all.ravel(), (Q_i.ravel(), Q_j.ravel())), 
                       shape=(nel * np_edof, sdof))

    invM_i = np.tile(np.arange(0, nel * np_edof, dtype=np.int32), (np_edof, 1)).T
    base_sequence = np.tile(np.arange(np_edof), nel * np_edof)
    offsets = np.repeat(np.arange(nel) * np_edof, np_edof**2)
    column_indices = base_sequence + offsets
    invM_all = csr_matrix((invM_all.ravel(), (invM_i.ravel(), column_indices.ravel())), 
                          shape=(nel * np_edof, nel * np_edof))

    # Compute pressure from incompressibility constraint
    Pressure = PF * invM_all @ (Q_all @ Vel)

    return Solution(Vel=Vel, Pressure=Pressure)
