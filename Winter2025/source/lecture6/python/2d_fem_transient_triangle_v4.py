# Script that calculated transient diffusion on unstructured triangle mesh
# V4: VECTORIZED BATCH PROCESSING - processes all elements simultaneously
# ============================================================================
# PERFORMANCE CONFIGURATION
# ============================================================================
USE_NUMBA = True

from dataclasses import dataclass
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from shapes_tri import shapes_tri
import triangle as tr
import meshio
import time

# Try to import numba
try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    if USE_NUMBA:
        print("Warning: Numba not available, falling back to pure Python")
    USE_NUMBA = False


@dataclass
class Mesh:
    """Container for mesh data."""
    GCOORD: np.ndarray  # Node coordinates, shape (nnod, 2)
    EL2NOD: np.ndarray  # Element connectivity, shape (nel, nnodel)
    Phases: np.ndarray  # Phase/material ID per element, shape (nel,)
    nnod: int = 0   # Number of nodes
    nel: int = 0    # Number of elements
    nnodel: int = 0 # Number of nodes per element


@dataclass
class GeometryParams:
    """Container for geometry parameters."""
    x0: float
    y0: float
    lx: float
    ly: float
    n_incl: int
    radius: float


@dataclass
class MaterialParams:
    """Container for material properties."""
    rho: float   # Density
    cp: float    # Heat capacity
    k1: float    # Thermal conductivity (matrix)
    k2: float    # Thermal conductivity (inclusions)


@dataclass
class BoundaryConditions:
    """Container for boundary conditions."""
    Ttop: float  # Temperature at top boundary
    Tbot: float  # Temperature at bottom boundary

@dataclass
class Variables:
    """Container for variables."""
    T: np.ndarray  # solution vector, shape (nnod,)
    Kel: np.ndarray  # element-wise thermal conductivity, shape (nel,)

@dataclass
class TimeParams:
    """Container for time stepping parameters."""
    dt: float         # Time step
    nt: int           # Number of time steps
    output_freq: int  # Output frequency (write every N timesteps)


@dataclass
class IntegrationParams:
    """Container for numerical integration parameters."""
    nip: int             # Number of integration points
    gauss: np.ndarray    # Gauss integration point coordinates, shape (nip, 2)
    weights: np.ndarray  # Integration weights, shape (nip,)
    N_all: np.ndarray = None      # Pre-computed shape functions, shape (nip, nnodel)
    dNds_all: np.ndarray = None   # Pre-computed derivatives, shape (nip, 2, nnodel)


def precompute_shape_functions(integration: IntegrationParams) -> IntegrationParams:
    """
    Pre-compute shape functions at all integration points.

    For triangular elements with fixed Gauss integration points, the shape
    functions N and their local derivatives dNds are constant at each
    integration point. Pre-computing eliminates repeated function calls.

    Parameters
    ----------
    integration : IntegrationParams
        Integration parameters with gauss points defined

    Returns
    -------
    IntegrationParams
        Updated integration parameters with N_all and dNds_all filled
    """
    nip = integration.nip
    nnodel = 3  # Triangular elements have 3 nodes

    # Allocate arrays for pre-computed values
    N_all = np.zeros((nip, nnodel))
    dNds_all = np.zeros((nip, 2, nnodel))

    # Compute shape functions at each integration point
    for ip in range(nip):
        xi = integration.gauss[ip, 0]
        eta = integration.gauss[ip, 1]
        N, dNds = shapes_tri(xi, eta)
        N_all[ip, :] = N
        dNds_all[ip, :, :] = dNds

    # Store in integration params
    integration.N_all = N_all
    integration.dNds_all = dNds_all

    return integration


def make_mesh(geom: GeometryParams) -> Mesh:
    """
    Create an unstructured triangle mesh with circular inclusions.

    Parameters
    ----------
    geom : GeometryParams
        Geometry parameters including domain size and inclusion properties

    Returns
    -------
    Mesh
        Mesh object containing node coordinates, connectivity, and phase IDs
    """
    ## Create the triangle mesh
    vertices = []
    segments = []
    regions = []

    # make a box with given dims and place given attribute at its center
    def _make_box(x: float, y: float, w: float, h: float, attribute: float,
                  vertices: list, segments: list, regions: list) -> None:
        # we modify the incoming lists
        i = len(vertices)

        vertices.extend([[x,   y],
                        [x+w, y],
                        [x+w, y+h],
                        [x,   y+h]])

        segments.extend([(i+0, i+1),
                        (i+1, i+2),
                        (i+2, i+3),
                        (i+3, i+0)])

        regions.append([x+0.01*w, y+0.01*h, attribute, 0.0001])

    def _make_inclusion(center_x: float, center_y: float, radius: float,
                       points_inc: int, attribute: float,
                       vertices: list, segments: list, regions: list) -> None:
        theta = np.linspace(0, 2*np.pi, points_inc, endpoint=False)
        xx = np.cos(theta)
        yy = np.sin(theta)

        i = len(vertices)

        vertices.extend(np.array([center_x + radius*xx, center_y + radius*yy]).T)

        Tmp = np.array([np.arange(i, i+points_inc), np.arange(i+1, i+points_inc+1)]).T
        Tmp[-1, 1] = i
        segments.extend(Tmp)

        regions.append([center_x, center_y, attribute, 0.0001])

    # generate input
    _make_box(geom.x0, geom.y0, geom.lx, geom.ly, 1, vertices, segments, regions)

    _make_inclusion(-0.8, -0.3, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(-0.5, -0.75, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(-0.6, 0.5, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(-0.1, -0.3, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(0.1, 0, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(0.5, -0.2, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(0.6, .3, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(0.7, .8, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(0, .75, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(-0.5, .05, geom.radius, 20, 100, vertices, segments, regions)
    _make_inclusion(0.5, -.75, geom.radius, 20, 100, vertices, segments, regions)

    A = dict(vertices=vertices, segments=segments, regions=regions)
    B = tr.triangulate(A, 'pq33Aa')

    # extract mesh information
    GCOORD = B.get("vertices")
    EL2NOD = B.get("triangles")
    Phases = B.get("triangle_attributes")

    # reshape Phases to 1D array
    nel = EL2NOD.shape[0]
    Phases = np.reshape(Phases, nel)

    return Mesh(GCOORD=GCOORD, EL2NOD=EL2NOD, Phases=Phases, nnod=GCOORD.shape[0], nel=nel, nnodel=EL2NOD.shape[1])


# ============================================================================
# VECTORIZED BATCH ASSEMBLY FUNCTIONS
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
    # Extract components
    a = matrices[:, 0, 0]
    b = matrices[:, 0, 1]
    c = matrices[:, 1, 0]
    d = matrices[:, 1, 1]
    
    # Compute determinant
    det = a * d - b * c
    
    # Compute inverse
    inv = np.zeros_like(matrices)
    inv[:, 0, 0] = d / det
    inv[:, 0, 1] = -b / det
    inv[:, 1, 0] = -c / det
    inv[:, 1, 1] = a / det
    
    return inv, det


def assemble_system_vectorized(
    GCOORD: np.ndarray,
    EL2NOD: np.ndarray,
    Kel: np.ndarray,
    T: np.ndarray,
    N_all: np.ndarray,
    dNds_all: np.ndarray,
    weights: np.ndarray,
    rho: float,
    cp: float,
    dt: float
) -> tuple:
    """
    VECTORIZED element assembly - processes all elements simultaneously.
    
    KEY OPTIMIZATION: Eliminates Python loop over elements by using
    vectorized NumPy operations and einsum for tensor contractions.
    
    Works for DEFORMING meshes - computes geometry every call.
    
    Parameters
    ----------
    GCOORD : np.ndarray
        Node coordinates, shape (nnod, 2)
    EL2NOD : np.ndarray
        Element connectivity, shape (nel, nnodel)
    Kel : np.ndarray
        Element conductivities, shape (nel,)
    T : np.ndarray
        Temperature field, shape (nnod,)
    N_all : np.ndarray
        Pre-computed shape functions, shape (nip, nnodel)
    dNds_all : np.ndarray
        Pre-computed derivatives, shape (nip, 2, nnodel)
    weights : np.ndarray
        Integration weights, shape (nip,)
    rho : float
        Density
    cp : float
        Heat capacity
    dt : float
        Time step
    
    Returns
    -------
    I, J, K : np.ndarray
        Sparse matrix triplet format
    Rhs_all : np.ndarray
        Right-hand side vector
    """
    nel = EL2NOD.shape[0]
    nnodel = EL2NOD.shape[1]
    nnod = GCOORD.shape[0]
    nip = N_all.shape[0]
    
    # Extract ALL element coordinates at once: shape (nel, nnodel, 2)
    ECOORD_all = GCOORD[EL2NOD]
    
    # Extract temperature at all elements: shape (nel, nnodel)
    T_all = T[EL2NOD]
    
    # Storage for element matrices: shape (nel, nnodel, nnodel)
    Ael_all = np.zeros((nel, nnodel, nnodel))
    Rhs_el_all = np.zeros((nel, nnodel))
    
    # Loop over integration points (small loop, only 3 iterations)
    for ip in range(nip):
        # Get pre-computed shape functions for this integration point
        N = N_all[ip, :]              # shape (nnodel,)
        dNds = dNds_all[ip, :, :]     # shape (2, nnodel)
        
        # Compute Jacobians for ALL elements at once using einsum
        # For each element: Jac = dNds @ ECOORD^T
        # dNds: (2, 3), ECOORD_all: (nel, 3, 2) -> Jac: (nel, 2, 2)
        # We want: Jac[n,i,j] = sum_k dNds[i,k] * ECOORD_all[n,k,j]
        Jac_all = np.einsum('ik,nkj->nij', dNds, ECOORD_all)
        
        # Vectorized inversion of all 2x2 Jacobians
        invJ_all, detJ_all = batch_invert_2x2(Jac_all)
        
        # Compute global derivatives for ALL elements: dNdx = invJ @ dNds
        # shape: (nel, 2, nnodel)
        dNdx_all = np.einsum('nij,jk->nik', invJ_all, dNds)
        
        # Mass matrix (same for all elements): Me = N ⊗ N
        # shape: (nnodel, nnodel)
        Me = np.outer(N, N)
        
        # Stiffness matrices for ALL elements: Ke = dNdx.T @ dNdx
        # dNdx_all: (nel, 2, nnodel), we want (nel, nnodel, nnodel)
        # Ke[n,i,j] = sum_k dNdx[n,k,i] * dNdx[n,k,j]
        Ke_all = np.einsum('nki,nkj->nij', dNdx_all, dNdx_all)
        
        # Integration weight times determinant for all elements
        # shape: (nel,)
        weight_detJ = weights[ip] * detJ_all
        
        # Assemble element matrices for ALL elements
        # Broadcasting: (nel, 1, 1) * (nnodel, nnodel) + (nel, 1, 1) * (nel, nnodel, nnodel)
        # shape: (nel, nnodel, nnodel)
        Ael_all += (rho * cp * Me[np.newaxis, :, :] + 
                    dt * Kel[:, np.newaxis, np.newaxis] * Ke_all) * weight_detJ[:, np.newaxis, np.newaxis]
        
        # Assemble RHS for ALL elements: Me @ T_el
        # shape: (nel, nnodel)
        Rhs_el_all += rho * cp * (Me @ T_all.T).T * weight_detJ[:, np.newaxis]
    
    # Build sparse matrix storage in triplet format
    # Create index arrays for sparse matrix
    I = np.zeros((nel, nnodel * nnodel), dtype=np.int32)
    J = np.zeros((nel, nnodel * nnodel), dtype=np.int32)
    K = np.zeros((nel, nnodel * nnodel))
    
    # Vectorized indexing
    for i in range(nnodel):
        for j in range(nnodel):
            idx = i * nnodel + j
            I[:, idx] = EL2NOD[:, i]
            J[:, idx] = EL2NOD[:, j]
            K[:, idx] = Ael_all[:, i, j]
    
    # Assemble global RHS by scattering element contributions
    Rhs_all = np.zeros(nnod)
    for i in range(nnodel):
        np.add.at(Rhs_all, EL2NOD[:, i], Rhs_el_all[:, i])
    
    return I, J, K, Rhs_all


def compute_heat_flux_vectorized(
    GCOORD: np.ndarray,
    EL2NOD: np.ndarray,
    Kel: np.ndarray,
    T: np.ndarray,
    N_centroid: np.ndarray,
    dNds_centroid: np.ndarray
) -> tuple:
    """
    VECTORIZED heat flux computation - processes all elements simultaneously.
    
    Parameters
    ----------
    GCOORD : np.ndarray
        Node coordinates, shape (nnod, 2)
    EL2NOD : np.ndarray
        Element connectivity, shape (nel, nnodel)
    Kel : np.ndarray
        Element conductivities, shape (nel,)
    T : np.ndarray
        Temperature field, shape (nnod,)
    N_centroid : np.ndarray
        Shape functions at centroid, shape (nnodel,)
    dNds_centroid : np.ndarray
        Derivatives at centroid, shape (2, nnodel)
    
    Returns
    -------
    Q_x, Q_y : np.ndarray
        Heat flux components, shape (nel,)
    """
    nel = EL2NOD.shape[0]
    
    # Extract ALL element coordinates: shape (nel, nnodel, 2)
    ECOORD_all = GCOORD[EL2NOD]
    
    # Extract temperature at all elements: shape (nel, nnodel)
    T_all = T[EL2NOD]
    
    # Compute Jacobians for ALL elements: shape (nel, 2, 2)
    # For each element: Jac = dNds @ ECOORD^T
    # dNds: (2, 3), ECOORD_all: (nel, 3, 2) -> Jac: (nel, 2, 2)
    Jac_all = np.einsum('ik,nkj->nij', dNds_centroid, ECOORD_all)
    
    # Vectorized inversion
    invJ_all, _ = batch_invert_2x2(Jac_all)
    
    # Global derivatives for ALL elements: shape (nel, 2, nnodel)
    dNdx_all = np.einsum('nij,jk->nik', invJ_all, dNds_centroid)
    
    # Heat flux for ALL elements: Q = -k * dNdx @ T
    # shape: (nel,)
    Q_x = -Kel * np.einsum('ni,ni->n', dNdx_all[:, 0, :], T_all)
    Q_y = -Kel * np.einsum('ni,ni->n', dNdx_all[:, 1, :], T_all)
    
    return Q_x, Q_y


# ============================================================================
# NUMBA-ACCELERATED VERSIONS (for comparison/fallback)
# ============================================================================

def assemble_system_numba(
    GCOORD: np.ndarray,
    EL2NOD: np.ndarray,
    Kel: np.ndarray,
    T: np.ndarray,
    N_all: np.ndarray,
    dNds_all: np.ndarray,
    weights: np.ndarray,
    rho: float,
    cp: float,
    dt: float
) -> tuple:
    """
    Element-by-element assembly with numba JIT compilation.
    Kept as fallback/comparison to vectorized version.
    """
    nel = EL2NOD.shape[0]
    nnodel = EL2NOD.shape[1]
    nnod = GCOORD.shape[0]
    nip = N_all.shape[0]

    # Storage
    Rhs_all = np.zeros(nnod)
    I = np.zeros((nel, nnodel * nnodel), dtype=np.int32)
    J = np.zeros((nel, nnodel * nnodel), dtype=np.int32)
    K = np.zeros((nel, nnodel * nnodel))

    for iel in range(nel):
        # Get element coordinates
        ECOORD = GCOORD[EL2NOD[iel, :], :]
        Ael = np.zeros((nnodel, nnodel))
        Rhs_el = np.zeros(nnodel)

        for ip in range(nip):
            # Get pre-computed shape functions
            N = N_all[ip, :]
            dNds = dNds_all[ip, :, :]

            # Jacobian, inverse, and determinant
            Jac = dNds @ ECOORD
            invJ = np.linalg.inv(Jac)
            detJ = np.linalg.det(Jac)

            # Global derivatives
            dNdx = invJ @ dNds

            # Element matrices
            Me = np.outer(N, N)
            Ke_local = dNdx.T @ dNdx

            # Assemble element matrix
            Ael += (rho * cp * Me + dt * Kel[iel] * Ke_local) * detJ * weights[ip]

            # Assemble right-hand side
            T_el = T[EL2NOD[iel, :]]
            Rhs_el += rho * cp * (Me @ T_el) * detJ * weights[ip]

        # Store coefficients for sparse matrix
        for i in range(nnodel):
            for j in range(nnodel):
                idx = i * nnodel + j
                I[iel, idx] = EL2NOD[iel, i]
                J[iel, idx] = EL2NOD[iel, j]
                K[iel, idx] = Ael[i, j]

        # Add to global RHS
        for i in range(nnodel):
            Rhs_all[EL2NOD[iel, i]] += Rhs_el[i]

    return I, J, K, Rhs_all


# ============================================================================
# CHOOSE ASSEMBLY FUNCTION BASED ON NUMBA AVAILABILITY
# ============================================================================

if USE_NUMBA and NUMBA_AVAILABLE:
    print("Numba available - compiling element-by-element assembly...")
    assemble_system_numba_jit = njit(cache=True)(assemble_system_numba)
    
    # Use vectorized as default (usually faster for large meshes)
    assemble_system_optimized = assemble_system_vectorized
    compute_heat_flux_optimized = compute_heat_flux_vectorized
    mode_name = "VECTORIZED"
else:
    # Use vectorized version (works without numba)
    assemble_system_optimized = assemble_system_vectorized
    compute_heat_flux_optimized = compute_heat_flux_vectorized
    mode_name = "VECTORIZED (no numba)"


def solve_2d_temperature_fem(
    mesh: Mesh,
    material: MaterialParams,
    bc: BoundaryConditions,
    time_params: TimeParams,
    integration: IntegrationParams,
    geom: GeometryParams,
    vars: Variables,
    writer: meshio.xdmf.TimeSeriesWriter,
    t: int
) -> Variables:
    """
    Solve 2D transient heat diffusion using FEM on triangular mesh.

    Parameters
    ----------
    mesh : Mesh
        Mesh object containing coordinates, connectivity, and phases
    material : MaterialParams
        Material properties (density, heat capacity, conductivities)
    bc : BoundaryConditions
        Boundary condition values (top and bottom temperatures)
    time_params : TimeParams
        Time stepping parameters
    integration : IntegrationParams
        Numerical integration parameters (Gauss points and weights)
    geom : GeometryParams
        Geometry parameters for boundary identification
    vars : Variables
        Variables container with temperature and conductivity
    writer : meshio.xdmf.TimeSeriesWriter
        Output file writer for visualization
    t : int
        Current time step index

    Returns
    -------
    vars : Variables
        Updated variables container with the temperature field
    """
    T = vars.T
    Kel = vars.Kel

    t_start_total = time.perf_counter()

    nnodel = mesh.EL2NOD.shape[1]
    nel = mesh.EL2NOD.shape[0]
    nnod = mesh.GCOORD.shape[0]

    # Element assembly - VECTORIZED!
    t_start_assembly = time.perf_counter()
    I, J, K, Rhs_all = assemble_system_optimized(
        mesh.GCOORD,
        mesh.EL2NOD,
        Kel,
        T,
        integration.N_all,
        integration.dNds_all,
        integration.weights,
        material.rho,
        material.cp,
        time_params.dt
    )
    t_assembly = time.perf_counter() - t_start_assembly

    # Create sparse matrix
    t_start_sparse = time.perf_counter()
    A_all = csr_matrix((K.reshape(nel*nnodel*nnodel), (I.reshape(nel*nnodel*nnodel), J.reshape(nel*nnodel*nnodel))), shape=(nnod, nnod))
    t_sparse = time.perf_counter() - t_start_sparse

    # Apply boundary conditions
    t_start_bc = time.perf_counter()
    # indices and values at top and bottom
    tol = 1e-3
    i_bot = np.where(abs(mesh.GCOORD[:, 1] - geom.y0) < tol)[0]
    i_top = np.where(abs(mesh.GCOORD[:, 1] - (geom.y0+geom.lx)) < tol)[0]

    Ind_bc = np.concatenate((i_bot, i_top))
    Val_bc = np.concatenate((np.ones(i_bot.shape)*bc.Tbot, np.ones(i_top.shape)*bc.Ttop))

    # boundary conditions that keep matrix symmetry
    Free = np.arange(0, nnod)
    Free = np.delete(Free, Ind_bc)
    TMP = A_all[:, Ind_bc]
    Rhs_all = Rhs_all - TMP.dot(Val_bc)
    t_bc = time.perf_counter() - t_start_bc

    # Solve linear system
    t_start_solve = time.perf_counter()
    T[Free] = spsolve(A_all[np.ix_(Free, Free)], Rhs_all[Free])
    T[Ind_bc] = Val_bc
    t_solve = time.perf_counter() - t_start_solve

    # Check if we should write output at this timestep
    write_output = (t % time_params.output_freq == 0) or (t == time_params.nt - 1)

    # Postprocessing - heat flow (ON-DEMAND: only when writing output)
    if write_output:
        t_start_postproc = time.perf_counter()

        # Pre-compute shape functions at element centroid
        N_centroid, dNds_centroid = shapes_tri(1/3, 1/3)

        # Compute heat flux - VECTORIZED!
        Q_x, Q_y = compute_heat_flux_optimized(
            mesh.GCOORD,
            mesh.EL2NOD,
            Kel,
            T,
            N_centroid,
            dNds_centroid
        )

        # Compute element centroids - vectorized
        ECOORD_all = mesh.GCOORD[mesh.EL2NOD]
        Ec_x = np.mean(ECOORD_all[:, :, 0], axis=1)
        Ec_y = np.mean(ECOORD_all[:, :, 1], axis=1)

        t_postproc = time.perf_counter() - t_start_postproc

        # Save data
        t_start_io = time.perf_counter()
        # cell data
        U = np.hstack((Q_x.reshape(-1, 1), Q_y.reshape(-1, 1)))
        U = np.hstack((U, U[:, 0].reshape(-1, 1)*0))

        # save data
        writer.write_data(t, point_data={"T": T}, cell_data={"U": [U], "K": [Kel]})
        t_io = time.perf_counter() - t_start_io
    else:
        # No output this timestep - skip postprocessing and I/O
        t_postproc = 0.0
        t_io = 0.0

    t_total = time.perf_counter() - t_start_total

    # Print timing information every 10 steps or on first step
    if t == 0 or (t + 1) % 10 == 0:
        output_str = " [OUTPUT WRITTEN]" if write_output else ""
        print(f"\n{'='*70}")
        print(f"Time step {t+1}/{time_params.nt} - Performance breakdown (Mode: {mode_name}){output_str}:")
        print(f"{'='*70}")
        print(f"  Element assembly:    {t_assembly*1000:8.2f} ms ({t_assembly/t_total*100:5.1f}%)")
        print(f"  Sparse matrix:       {t_sparse*1000:8.2f} ms ({t_sparse/t_total*100:5.1f}%)")
        print(f"  Boundary conditions: {t_bc*1000:8.2f} ms ({t_bc/t_total*100:5.1f}%)")
        print(f"  Linear solve:        {t_solve*1000:8.2f} ms ({t_solve/t_total*100:5.1f}%)")
        if write_output:
            print(f"  Postprocessing:      {t_postproc*1000:8.2f} ms ({t_postproc/t_total*100:5.1f}%)")
            print(f"  File I/O:            {t_io*1000:8.2f} ms ({t_io/t_total*100:5.1f}%)")
        else:
            print(f"  Postprocessing:      SKIPPED (on-demand)")
            print(f"  File I/O:            SKIPPED (on-demand)")
        print(f"  {'─'*70}")
        print(f"  Total time:          {t_total*1000:8.2f} ms")
        print(f"  Mesh info: {nel} elements, {nnod} nodes, {integration.nip} integration points/element")
        print(f"{'='*70}")

    vars.T = T
    return vars


def main() -> None:
    """
    Main driver for 2D transient heat diffusion FEM simulation.
    """
    print("\n" + "="*70)
    print("2D TRANSIENT HEAT DIFFUSION - FEM SOLVER V4")
    print("VECTORIZED BATCH PROCESSING - All elements computed simultaneously")
    print("="*70)

    # Geometry parameters
    geom = GeometryParams(
        x0=-1.0,
        y0=-1.0,
        lx=2.0,
        ly=2.0,
        n_incl=5,
        radius=0.15
    )

    # Time parameters
    time_params = TimeParams(
        dt=0.015,       # arbitrary time stepping
        nt=80,          # number of timesteps
        output_freq=10  # write output every 10 timesteps
    )

    # Material parameters
    material = MaterialParams(
        rho=1.0,
        cp=1.0,
        k1=1.0,
        k2=0.01
    )

    # Boundary conditions
    bc = BoundaryConditions(
        Ttop=0.0,
        Tbot=1.0
    )

    # Gauss integration points for triangles
    integration = IntegrationParams(
        nip=3,
        gauss=np.array([[1/6, 2/3, 1/6], [1/6, 1/6, 2/3]]).T.copy(),
        weights=np.array([1/6, 1/6, 1/6])
    )

    # Pre-compute shape functions at integration points
    print("\nPre-computing shape functions at integration points...")
    integration = precompute_shape_functions(integration)
    print(f"  - Shape functions pre-computed for {integration.nip} integration points")

    # variables
    vars = Variables(
        T=None,
        Kel=None
    )

    # Create mesh
    print("\nGenerating mesh...")
    t_start_mesh = time.perf_counter()
    mesh = make_mesh(geom)
    t_mesh = time.perf_counter() - t_start_mesh

    # Extract mesh information
    nel = mesh.EL2NOD.shape[0]
    nnod = mesh.GCOORD.shape[0]

    print(f"Mesh generated in {t_mesh*1000:.2f} ms")
    print(f"  - Number of elements: {nel}")
    print(f"  - Number of nodes: {nnod}")
    print(f"  - Nodes per element: {mesh.nnodel}")

    # Setup element-wise thermal conductivity
    n_matrix = np.sum(mesh.Phases == 1)
    n_inclusions = np.sum(mesh.Phases == 100)
    print(f"\nMaterial distribution:")
    print(f"  - Matrix elements (k={material.k1}): {n_matrix}")
    print(f"  - Inclusion elements (k={material.k2}): {n_inclusions}")

    Kel = np.ones(nel) * material.k1
    Kel[np.where(mesh.Phases == 100)] = material.k2
    vars.Kel = Kel
    vars.T = np.zeros(nnod) # initial T

    # Setup output writing
    print(f"\nSetting up output files...")
    points = np.hstack((mesh.GCOORD, mesh.GCOORD[:, 0].reshape(-1, 1)*0))  # must have 3 components (x,y,z)
    cells = [("triangle", mesh.EL2NOD)]
    writer = meshio.xdmf.TimeSeriesWriter('transient_v4.xmf')
    writer.__enter__()  # have to add this: import hdf5 and open file ...
    writer.write_points_cells(points, cells)

    # Time loop
    num_outputs = (time_params.nt - 1) // time_params.output_freq + 1  # +1 for final step
    print(f"\nStarting time integration ({time_params.nt} steps, dt={time_params.dt})...")
    print(f"  - Output frequency: every {time_params.output_freq} steps ({num_outputs} outputs total)")
    print(f"  - Assembly mode: {mode_name}")
    t_start_timeloop = time.perf_counter()

    for t in range(0, time_params.nt):
        vars = solve_2d_temperature_fem(
            mesh=mesh,
            material=material,
            bc=bc,
            time_params=time_params,
            integration=integration,
            geom=geom,
            vars=vars,
            writer=writer,
            t=t
        )

    t_timeloop = time.perf_counter() - t_start_timeloop

    writer.__exit__()  # close file

    # Final summary
    print(f"\n{'='*70}")
    print("SIMULATION COMPLETE - PERFORMANCE SUMMARY")
    print(f"{'='*70}")
    print(f"Assembly mode:         {mode_name}")
    print(f"Mesh generation:       {t_mesh*1000:10.2f} ms")
    print(f"Time integration:      {t_timeloop:10.2f} s  ({time_params.nt} steps)")
    print(f"Average per step:      {t_timeloop/time_params.nt*1000:10.2f} ms")
    print(f"Total simulation time: {(t_mesh + t_timeloop):10.2f} s")
    print(f"{'='*70}")
    print(f"Output files: transient_v4.xmf, transient_v4.h5")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
