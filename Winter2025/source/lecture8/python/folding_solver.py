"""
2D Mechanical FEM solver for Stokes flow with viscous folding support.

Based on mechanical2d_blocked.py (blocked vectorized processing).
Adds time-stepping infrastructure for Lagrangian advection of a folding layer.
"""

# ============================================================================
# PERFORMANCE CONFIGURATION
# ============================================================================
BLOCK_SIZE = 10000  # Process elements in blocks for better cache utilization

from dataclasses import dataclass
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from shp_deriv_triangle import shp_deriv_triangle
from ip_triangle import ip_triangle
import time


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


@dataclass
class TimeParams:
    """Container for time-stepping parameters."""
    dt: float           # Time step size
    nt: int             # Number of time steps
    output_freq: int    # Output every N steps


@dataclass
class FoldingParams:
    """Container for folding-specific parameters."""
    Lx: float           # Domain width
    Ly: float           # Domain height
    h: float            # Layer thickness
    A: float            # Perturbation amplitude
    wavelength: float   # Perturbation wavelength
    mu_matrix: float    # Matrix viscosity
    mu_layer: float     # Layer viscosity
    eps_bg: float       # Background strain rate
    el_size_layer: float  # Element size in/near layer
    el_size_matrix: float # Element size in matrix
    n_interface_pts: int = 100  # Points along each interface


# ============================================================================
# BLOCKED VECTORIZED HELPER FUNCTIONS
# ============================================================================

def batch_invert_2x2(matrices: np.ndarray) -> tuple:
    """
    Vectorized inversion of 2x2 matrices.

    Parameters
    ----------
    matrices : np.ndarray
        Array of 2x2 matrices, shape (n, 2, 2)

    Returns
    -------
    inv_matrices : np.ndarray
        Inverted matrices, shape (n, 2, 2)
    determinants : np.ndarray
        Determinants, shape (n,)
    """
    a = matrices[:, 0, 0]
    b = matrices[:, 0, 1]
    c = matrices[:, 1, 0]
    d = matrices[:, 1, 1]

    det = a * d - b * c

    inv = np.zeros_like(matrices)
    inv[:, 0, 0] = d / det
    inv[:, 0, 1] = -b / det
    inv[:, 1, 0] = -c / det
    inv[:, 1, 1] = a / det

    return inv, det


# ============================================================================
# LAGRANGIAN ADVECTION FUNCTIONS
# ============================================================================

def advect_mesh(mesh: Mesh, velocity: np.ndarray, dt: float) -> None:
    """
    Advect mesh nodes with the velocity field (Lagrangian update).

    Moves all real nodes (corner + edge) by velocity*dt, then recomputes
    center node positions as the mean of the three corner nodes.

    Parameters
    ----------
    mesh : Mesh
        Mesh object (modified in-place)
    velocity : np.ndarray
        Velocity field, shape (2*nnod,) interleaved [vx0, vy0, vx1, vy1, ...]
    dt : float
        Time step size
    """
    Vx = velocity[0::2]
    Vy = velocity[1::2]

    # Number of real nodes (corner + edge midpoints, excludes center nodes)
    n_real = mesh.nnod - mesh.nel

    # Move real nodes
    mesh.GCOORD[:n_real, 0] += Vx[:n_real] * dt
    mesh.GCOORD[:n_real, 1] += Vy[:n_real] * dt

    # Recompute center node positions (node 7 = mean of corner nodes 0-2)
    mesh.GCOORD[n_real:, :] = np.mean(
        mesh.GCOORD[mesh.EL2NOD[:, 0:3]], axis=1
    )


def update_boundary_conditions(mesh: Mesh, eps_bg: float) -> BoundaryConditions:
    """
    Recompute pure shear boundary conditions on the (deformed) mesh.

    Pure shear: vx = -eps_bg * x, vy = +eps_bg * y
    Applied on all boundary nodes (markers 101-104).

    Parameters
    ----------
    mesh : Mesh
        Current mesh (possibly deformed)
    eps_bg : float
        Background strain rate

    Returns
    -------
    BoundaryConditions
        Updated boundary conditions
    """
    bc_nodes = np.where(np.isin(mesh.Node_ids, [101, 102, 103, 104]))[0]

    Bc_val_x = -eps_bg * mesh.GCOORD[bc_nodes, 0]  # horizontal compression
    Bc_val_y =  eps_bg * mesh.GCOORD[bc_nodes, 1]  # vertical extension

    Bc_val = np.hstack((Bc_val_x, Bc_val_y))
    Bc_ind = np.hstack((2 * bc_nodes, 2 * bc_nodes + 1))

    return BoundaryConditions(Bc_ind=Bc_ind, Bc_val=Bc_val)


# ============================================================================
# BLOCKED VECTORIZED SOLVER
# ============================================================================

def solve_mechanical_2d(
    mesh: Mesh,
    material: MaterialParams,
    bc: BoundaryConditions,
    solver: SolverParams,
    block_size: int = BLOCK_SIZE
) -> Solution:
    """
    Solve 2D incompressible Stokes flow using FEM with discontinuous pressure.

    BLOCKED VECTORIZED assembly: processes elements in blocks for optimal
    cache usage while maintaining full NumPy vectorization within each block.

    Uses triangular elements with quadratic velocity and linear
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
    block_size : int
        Number of elements per block (default: 10000 for optimal cache usage)

    Returns
    -------
    Solution
        Solution object containing velocity and pressure fields
    """
    t_start_total = time.perf_counter()

    # Extract mesh dimensions
    nel = mesh.nel
    nnod = mesh.nnod
    nnodel = mesh.nnodel

    # CONSTANTS
    ndim = 2
    nedof = nnodel * ndim      # DOFs per element (7*2 = 14)
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
    N_all = Shapes['shape_functions']       # (nip, nnodel)
    dNds_all = Shapes['shape_func_deriv']   # (nip, 2, nnodel)

    # BUILD ELEMENT-TO-DOF MAPPING
    EL2DOF = np.zeros((nel, nedof), dtype=np.int32)
    EL2DOF[:, 0::ndim] = ndim * mesh.EL2NOD
    EL2DOF[:, 1::ndim] = ndim * mesh.EL2NOD + 1

    # ALLOCATE GLOBAL STORAGE
    A_all = np.zeros((nel, nedof * nedof))
    Q_all = np.zeros((nel, nedof * np_edof))
    invM_all = np.zeros((nel, np_edof * np_edof))
    Rhs_all = np.zeros(sdof)

    # BUILD SPARSE MATRIX INDICES (vectorized, outside block loop)
    I_all = np.repeat(EL2DOF, nedof, axis=1)    # each DOF repeated nedof times
    J_all = np.tile(EL2DOF, (1, nedof))          # DOF array repeated nedof times

    # Calculate number of blocks
    n_blocks = (nel + block_size - 1) // block_size

    # ================================================================
    # BLOCKED ELEMENT ASSEMBLY
    # ================================================================
    t_start_assembly = time.perf_counter()

    for block_idx in range(n_blocks):
        # Determine block range
        start_el = block_idx * block_size
        end_el = min((block_idx + 1) * block_size, nel)
        n_el_block = end_el - start_el

        # Extract block data
        EL2NOD_block = mesh.EL2NOD[start_el:end_el]
        EL2DOF_block = EL2DOF[start_el:end_el]
        ECOORD_block = mesh.GCOORD[EL2NOD_block]                   # (n_el, nnodel, 2)
        EMu_block = material.Mu[mesh.Phases[start_el:end_el]]      # (n_el,)
        ERho_block = material.Rho[mesh.Phases[start_el:end_el]]    # (n_el,)

        # Initialize block element matrices
        A_block = np.zeros((n_el_block, nedof, nedof))
        Q_block = np.zeros((n_el_block, np_edof, nedof))
        M_block = np.zeros((n_el_block, np_edof, np_edof))
        Rhs_block = np.zeros((n_el_block, nedof))

        # Pressure interpolation matrix P (constant per element)
        # P = [[1, 1, 1], [x0, x1, x2], [y0, y1, y2]]
        P_block = np.ones((n_el_block, np_edof, np_edof))
        P_block[:, 1:3, :] = ECOORD_block[:, :3, :].transpose(0, 2, 1)

        # Precompute P inverse (reused for each integration point)
        P_inv_block = np.linalg.inv(P_block)    # (n_el, 3, 3)

        # INTEGRATION POINT LOOP
        for ip in range(solver.nip):
            # Load pre-computed shape functions for this integration point
            Ni = N_all[ip, :]              # (nnodel,)
            dNdui = dNds_all[ip, :, :]     # (2, nnodel)

            # Compute pressure basis Pb at this integration point
            Pb_block = np.ones((n_el_block, np_edof))
            ip_coords = np.einsum('k,nkj->nj', Ni, ECOORD_block)   # (n_el, 2)
            Pb_block[:, 1] = ip_coords[:, 0]
            Pb_block[:, 2] = ip_coords[:, 1]

            # Pressure shape functions: Pi = P^{-1} @ Pb
            Pi_block = np.einsum('nij,nj->ni', P_inv_block, Pb_block)   # (n_el, 3)

            # Compute Jacobians for ALL elements in block
            Jac_block = np.einsum('ik,nkj->nij', dNdui, ECOORD_block)  # (n_el, 2, 2)

            # Vectorized inversion of all 2x2 Jacobians
            invJ_block, detJ_block = batch_invert_2x2(Jac_block)

            # Global derivatives of shape functions
            dNdX_block = np.einsum('nij,jk->nik', invJ_block, dNdui)   # (n_el, 2, nnodel)

            # Integration weight
            weight_block = IP_w[ip] * detJ_block    # (n_el,)

            # Build strain-displacement matrix B for ALL elements in block
            B_block = np.zeros((n_el_block, 3, nedof))
            B_block[:, 0, 0::2] = dNdX_block[:, 0, :]  # dNdx for x-displacement
            B_block[:, 1, 1::2] = dNdX_block[:, 1, :]  # dNdy for y-displacement
            B_block[:, 2, 0::2] = dNdX_block[:, 1, :]  # dNdy for x-displacement (shear)
            B_block[:, 2, 1::2] = dNdX_block[:, 0, :]  # dNdx for y-displacement (shear)

            # Volumetric strain: Bvol = dNdX^T, raveled as [dNdx_0, dNdy_0, ...]
            Bvol_ravel_block = dNdX_block.transpose(0, 2, 1).reshape(n_el_block, -1)  # (n_el, nedof)

            # B^T @ DEV @ B for ALL elements in block
            DEV_B = np.einsum('ij,njk->nik', DEV, B_block)             # (n_el, 3, nedof)
            BtDB_block = np.einsum('nji,njk->nik', B_block, DEV_B)     # (n_el, nedof, nedof)

            # Accumulate element stiffness matrices
            A_block += (weight_block * EMu_block)[:, np.newaxis, np.newaxis] * BtDB_block

            # Accumulate pressure-velocity coupling Q
            Q_block -= weight_block[:, np.newaxis, np.newaxis] * np.einsum('ni,nj->nij', Pi_block, Bvol_ravel_block)

            # Accumulate pressure mass matrix M
            M_block += weight_block[:, np.newaxis, np.newaxis] * np.einsum('ni,nj->nij', Pi_block, Pi_block)

            # Body force (same pattern for all elements, only weight*rho varies)
            body_force = np.outer(material.G, Ni).T.ravel()     # (nedof,)
            Rhs_block += (weight_block * ERho_block)[:, np.newaxis] * body_force[np.newaxis, :]

        # STATIC CONDENSATION of pressure DOFs
        invM_block = np.linalg.inv(M_block)                         # (n_el, 3, 3)
        invM_Q = np.einsum('nij,njk->nik', invM_block, Q_block)    # (n_el, 3, nedof)
        QtinvMQ = np.einsum('nji,njk->nik', Q_block, invM_Q)       # (n_el, nedof, nedof)
        A_block += PF * QtinvMQ

        # Store element contributions
        A_all[start_el:end_el, :] = A_block.reshape(n_el_block, -1)
        Q_all[start_el:end_el, :] = Q_block.reshape(n_el_block, -1)
        invM_all[start_el:end_el, :] = invM_block.reshape(n_el_block, -1)

        # Scatter RHS contributions to global vector
        for i in range(nedof):
            np.add.at(Rhs_all, EL2DOF_block[:, i], Rhs_block[:, i])

    t_assembly = time.perf_counter() - t_start_assembly

    # ================================================================
    # ASSEMBLE GLOBAL SPARSE MATRIX
    # ================================================================
    t_start_sparse = time.perf_counter()
    A_sparse = csr_matrix((A_all.ravel(), (I_all.ravel(), J_all.ravel())),
                          shape=(sdof, sdof))
    t_sparse = time.perf_counter() - t_start_sparse

    # ================================================================
    # APPLY BOUNDARY CONDITIONS (preserving symmetry)
    # ================================================================
    t_start_bc = time.perf_counter()
    Free = np.arange(0, sdof)
    Free = np.delete(Free, bc.Bc_ind)
    TMP = A_sparse[:, bc.Bc_ind]
    Rhs_all = Rhs_all - TMP.dot(bc.Bc_val)
    t_bc = time.perf_counter() - t_start_bc

    # ================================================================
    # SOLVE LINEAR SYSTEM
    # ================================================================
    t_start_solve = time.perf_counter()
    Vel = np.zeros(sdof)
    Vel[Free] = spsolve(A_sparse[np.ix_(Free, Free)], Rhs_all[Free])
    Vel[bc.Bc_ind] = bc.Bc_val
    t_solve = time.perf_counter() - t_start_solve

    # ================================================================
    # RECOVER PRESSURE FIELD
    # ================================================================
    t_start_pressure = time.perf_counter()

    # Build global Q matrix for discontinuous pressure
    Q_i = np.tile(np.arange(0, nel * np_edof, dtype=np.int32), (nedof, 1)).T
    Q_j = np.tile(EL2DOF, (1, np_edof))
    Q_sparse = csr_matrix((Q_all.ravel(), (Q_i.ravel(), Q_j.ravel())),
                          shape=(nel * np_edof, sdof))

    # Build global invM matrix (block-diagonal)
    invM_i = np.tile(np.arange(0, nel * np_edof, dtype=np.int32), (np_edof, 1)).T
    base_sequence = np.tile(np.arange(np_edof), nel * np_edof)
    offsets = np.repeat(np.arange(nel) * np_edof, np_edof**2)
    column_indices = base_sequence + offsets
    invM_sparse = csr_matrix((invM_all.ravel(), (invM_i.ravel(), column_indices.ravel())),
                             shape=(nel * np_edof, nel * np_edof))

    # Compute pressure from incompressibility constraint
    Pressure = PF * invM_sparse @ (Q_sparse @ Vel)

    t_pressure = time.perf_counter() - t_start_pressure
    t_total = time.perf_counter() - t_start_total

    # ================================================================
    # PERFORMANCE REPORT
    # ================================================================
    print(f"\n{'='*70}")
    print(f"SOLVER PERFORMANCE - BLOCKED VECTORIZED (block_size={block_size})")
    print(f"{'='*70}")
    print(f"  Element assembly:    {t_assembly*1000:8.2f} ms ({t_assembly/t_total*100:5.1f}%)")
    print(f"  Sparse matrix:       {t_sparse*1000:8.2f} ms ({t_sparse/t_total*100:5.1f}%)")
    print(f"  Boundary conditions: {t_bc*1000:8.2f} ms ({t_bc/t_total*100:5.1f}%)")
    print(f"  Linear solve:        {t_solve*1000:8.2f} ms ({t_solve/t_total*100:5.1f}%)")
    print(f"  Pressure recovery:   {t_pressure*1000:8.2f} ms ({t_pressure/t_total*100:5.1f}%)")
    print(f"  {'─'*68}")
    print(f"  Total time:          {t_total*1000:8.2f} ms")
    print(f"  Mesh: {nel} elements, {nnod} nodes, {n_blocks} blocks")
    print(f"{'='*70}")

    return Solution(Vel=Vel, Pressure=Pressure)
