# Script that calculates Turing pattern formation using double diffusion on unstructured triangle mesh
# Optimized version with Numba JIT compilation and structured function approach
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
from numpy.random import default_rng
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
    nnod: int = 0       # Number of nodes
    nel: int = 0        # Number of elements
    nnodel: int = 0     # Number of nodes per element
    sdof: int = 0       # System degrees of freedom (2 per node)


@dataclass
class GeometryParams:
    """Container for geometry parameters."""
    x0: float
    y0: float
    lx: float
    ly: float


@dataclass
class ModelParams:
    """Container for Turing pattern model parameters."""
    g_coeff: float   # Kinetic constant
    a_coeff: float   # Growth parameter A
    b_coeff: float   # Growth parameter B
    d_coeff: float   # Diffusivity B
    amp: float       # Random noise amplitude


@dataclass
class TimeParams:
    """Container for time stepping parameters."""
    dt: float   # Time step
    nt: int     # Number of time steps


@dataclass
class SolverParams:
    """Container for nonlinear solver parameters."""
    tol: float      # Convergence tolerance
    iter_max: int   # Maximum iterations


@dataclass
class IntegrationParams:
    """Container for numerical integration parameters."""
    nip: int             # Number of integration points
    gauss: np.ndarray    # Gauss integration point coordinates, shape (nip, 2)
    weights: np.ndarray  # Integration weights, shape (nip,)
    N_all: np.ndarray = None      # Pre-computed shape functions, shape (nip, nnodel)
    dNds_all: np.ndarray = None   # Pre-computed derivatives, shape (nip, 2, nnodel)


@dataclass
class Variables:
    """Container for solution variables."""
    A: np.ndarray  # Species A concentration, shape (nnod,)
    B: np.ndarray  # Species B concentration, shape (nnod,)


def precompute_shape_functions(integration: IntegrationParams) -> IntegrationParams:
    """
    Pre-compute shape functions at all integration points.

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
    Create an unstructured triangle mesh.

    Parameters
    ----------
    geom : GeometryParams
        Geometry parameters including domain size

    Returns
    -------
    Mesh
        Mesh object containing node coordinates, connectivity, and phase IDs
    """
    vertices = []
    segments = []
    regions = []

    def _make_box(x: float, y: float, w: float, h: float, attribute: float,
                  vertices: list, segments: list, regions: list) -> None:
        i = len(vertices)

        vertices.extend([[x,   y],
                        [x+w, y],
                        [x+w, y+h],
                        [x,   y+h]])

        segments.extend([(i+0, i+1),
                        (i+1, i+2),
                        (i+2, i+3),
                        (i+3, i+0)])

        regions.append([x+0.01*w, y+0.01*h, attribute, 0.005])

    # Generate input
    _make_box(geom.x0, geom.y0, geom.lx, geom.ly, 1, vertices, segments, regions)
    A = dict(vertices=vertices, segments=segments, regions=regions)
    B = tr.triangulate(A, 'pq33Aa')

    # Extract mesh information
    GCOORD = B.get("vertices")
    EL2NOD = B.get("triangles")
    Phases = B.get("triangle_attributes")

    nel = EL2NOD.shape[0]
    nnod = GCOORD.shape[0]
    nnodel = EL2NOD.shape[1]
    sdof = nnod * 2  # Two dof per node (A and B)
    Phases = np.reshape(Phases, nel)

    return Mesh(GCOORD=GCOORD, EL2NOD=EL2NOD, Phases=Phases,
                nnod=nnod, nel=nel, nnodel=nnodel, sdof=sdof)


def initialize_concentrations(mesh: Mesh, model: ModelParams) -> Variables:
    """
    Initialize species concentrations with random perturbations.

    Parameters
    ----------
    mesh : Mesh
        Mesh object
    model : ModelParams
        Model parameters

    Returns
    -------
    Variables
        Initial concentration fields
    """
    rng = default_rng()
    vals = rng.standard_normal(mesh.nnod)
    A = (model.a_coeff + model.b_coeff) + model.amp * vals

    vals = rng.standard_normal(mesh.nnod)
    B = model.b_coeff / (model.a_coeff + model.b_coeff)**2 + model.amp * vals

    return Variables(A=A, B=B)


# ============================================================================
# ASSEMBLY FUNCTIONS (Same code, optionally JIT-compiled by numba)
# ============================================================================

def assemble_system_optimized(
    GCOORD: np.ndarray,
    EL2NOD: np.ndarray,
    A: np.ndarray,
    B: np.ndarray,
    N_all: np.ndarray,
    dNds_all: np.ndarray,
    weights: np.ndarray,
    g_coeff: float,
    d_coeff: float,
    dt: float
) -> tuple:
    """
    Assemble stiffness matrices and RHS for coupled diffusion system.

    This function works as pure Python OR can be JIT-compiled by numba.
    The code is IDENTICAL in both cases - only execution speed differs.

    Parameters
    ----------
    GCOORD : np.ndarray
        Node coordinates, shape (nnod, 2)
    EL2NOD : np.ndarray
        Element connectivity, shape (nel, nnodel)
    A : np.ndarray
        Species A concentration, shape (nnod,)
    B : np.ndarray
        Species B concentration, shape (nnod,)
    N_all : np.ndarray
        Pre-computed shape functions, shape (nip, nnodel)
    dNds_all : np.ndarray
        Pre-computed derivatives, shape (nip, 2, nnodel)
    weights : np.ndarray
        Integration weights, shape (nip,)
    g_coeff : float
        Kinetic constant
    d_coeff : float
        Diffusivity for species B
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
    sdof = nnod * 2

    # Storage
    Rhs_all = np.zeros(sdof)
    I = np.zeros((nel, 2 * nnodel * nnodel))
    J = np.zeros((nel, 2 * nnodel * nnodel))
    K = np.zeros((nel, 2 * nnodel * nnodel))

    for iel in range(nel):
        # Get element coordinates
        ECOORD = GCOORD[EL2NOD[iel, :], :]
        Ael_A = np.zeros((nnodel, nnodel))
        Ael_B = np.zeros((nnodel, nnodel))
        RhsA_el = np.zeros(nnodel)
        RhsB_el = np.zeros(nnodel)

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

            # Mass and stiffness components
            Me = np.outer(N, N)
            Ke = dNdx.T @ dNdx

            # Element matrices for species A and B
            Ael_A += (Me * (1 + g_coeff * dt) + dt * Ke) * detJ * weights[ip]
            Ael_B += (Me + d_coeff * dt * Ke) * detJ * weights[ip]

            # RHS assembly (initial condition contribution)
            A_el = A[EL2NOD[iel, :]]
            B_el = B[EL2NOD[iel, :]]
            RhsA_el += (Me @ A_el) * detJ * weights[ip]
            RhsB_el += (Me @ B_el) * detJ * weights[ip]

        # Store coefficients for sparse matrix (interleaved DOFs)
        # A: even DOFs (0, 2, 4, ...), B: odd DOFs (1, 3, 5, ...)
        for i in range(nnodel):
            for j in range(nnodel):
                idx_A = i * nnodel + j
                idx_B = nnodel * nnodel + i * nnodel + j
                # A equations
                I[iel, idx_A] = 2 * EL2NOD[iel, i]
                J[iel, idx_A] = 2 * EL2NOD[iel, j]
                K[iel, idx_A] = Ael_A[i, j]
                # B equations
                I[iel, idx_B] = 2 * EL2NOD[iel, i] + 1
                J[iel, idx_B] = 2 * EL2NOD[iel, j] + 1
                K[iel, idx_B] = Ael_B[i, j]

        # Add to global RHS (interleaved)
        for i in range(nnodel):
            Rhs_all[2 * EL2NOD[iel, i]] += RhsA_el[i]
            Rhs_all[2 * EL2NOD[iel, i] + 1] += RhsB_el[i]

    return I, J, K, Rhs_all


def assemble_nonlinear_rhs_optimized(
    GCOORD: np.ndarray,
    EL2NOD: np.ndarray,
    A: np.ndarray,
    B: np.ndarray,
    N_all: np.ndarray,
    dNds_all: np.ndarray,
    weights: np.ndarray,
    g_coeff: float,
    a_coeff: float,
    b_coeff: float,
    dt: float
) -> np.ndarray:
    """
    Assemble nonlinear reaction terms for the RHS.

    This function works as pure Python OR can be JIT-compiled by numba.
    The code is IDENTICAL in both cases - only execution speed differs.

    Parameters
    ----------
    GCOORD : np.ndarray
        Node coordinates, shape (nnod, 2)
    EL2NOD : np.ndarray
        Element connectivity, shape (nel, nnodel)
    A : np.ndarray
        Species A concentration, shape (nnod,)
    B : np.ndarray
        Species B concentration, shape (nnod,)
    N_all : np.ndarray
        Pre-computed shape functions, shape (nip, nnodel)
    dNds_all : np.ndarray
        Pre-computed derivatives, shape (nip, 2, nnodel)
    weights : np.ndarray
        Integration weights, shape (nip,)
    g_coeff : float
        Kinetic constant
    a_coeff : float
        Growth parameter A
    b_coeff : float
        Growth parameter B
    dt : float
        Time step

    Returns
    -------
    F_all : np.ndarray
        Nonlinear force vector
    """
    nel = EL2NOD.shape[0]
    nnodel = EL2NOD.shape[1]
    nnod = GCOORD.shape[0]
    nip = N_all.shape[0]
    sdof = nnod * 2

    F_all = np.zeros(sdof)

    for iel in range(nel):
        ECOORD = GCOORD[EL2NOD[iel, :], :]
        FA_el = np.zeros(nnodel)
        FB_el = np.zeros(nnodel)

        for ip in range(nip):
            # Get pre-computed shape functions
            N = N_all[ip, :]
            dNds = dNds_all[ip, :, :]

            # Jacobian and determinant
            Jac = dNds @ ECOORD
            detJ = np.linalg.det(Jac)

            # Interpolate concentrations at integration point
            A_el = A[EL2NOD[iel, :]]
            B_el = B[EL2NOD[iel, :]]
            Ai = N @ A_el
            Bi = N @ B_el

            # Reaction terms: dA/dt = g*(a + B*A^2), dB/dt = g*(b - B*A^2)
            FA_el += N * dt * g_coeff * (a_coeff + Bi * Ai**2) * detJ * weights[ip]
            FB_el += N * dt * g_coeff * (b_coeff - Bi * Ai**2) * detJ * weights[ip]

        # Add to global force vector (interleaved)
        for i in range(nnodel):
            F_all[2 * EL2NOD[iel, i]] += FA_el[i]
            F_all[2 * EL2NOD[iel, i] + 1] += FB_el[i]

    return F_all


# ============================================================================
# APPLY NUMBA JIT COMPILATION (if enabled)
# ============================================================================
if USE_NUMBA and NUMBA_AVAILABLE:
    print("Applying numba JIT compilation to assembly functions...")
    assemble_system_optimized = njit(cache=True)(assemble_system_optimized)
    assemble_nonlinear_rhs_optimized = njit(cache=True)(assemble_nonlinear_rhs_optimized)


def solve_timestep(
    mesh: Mesh,
    model: ModelParams,
    time_params: TimeParams,
    solver: SolverParams,
    integration: IntegrationParams,
    vars: Variables,
    A_all: csr_matrix,
    Rhs_base: np.ndarray,
    t: int
) -> tuple:
    """
    Solve one time step with nonlinear iterations.

    Parameters
    ----------
    mesh : Mesh
        Mesh object
    model : ModelParams
        Model parameters
    time_params : TimeParams
        Time stepping parameters
    solver : SolverParams
        Solver parameters
    integration : IntegrationParams
        Integration parameters
    vars : Variables
        Current solution variables
    A_all : csr_matrix
        System matrix
    Rhs_base : np.ndarray
        Base RHS vector (before nonlinear terms)
    t : int
        Current time step

    Returns
    -------
    vars : Variables
        Updated solution variables
    iter : int
        Number of iterations
    error : float
        Final error
    """
    A = vars.A
    B = vars.B

    # Nonlinear iteration loop
    iter_count = 0
    error = 10.0
    Conc_tmp = np.ones(mesh.sdof) * 10

    t_start_iter = time.perf_counter()

    while error > solver.tol and iter_count < solver.iter_max:
        iter_count += 1

        # Assemble nonlinear RHS terms
        F_nl = assemble_nonlinear_rhs_optimized(
            mesh.GCOORD,
            mesh.EL2NOD,
            A,
            B,
            integration.N_all,
            integration.dNds_all,
            integration.weights,
            model.g_coeff,
            model.a_coeff,
            model.b_coeff,
            time_params.dt
        )

        # Complete RHS
        Rhs = Rhs_base + F_nl

        # Solve system
        Conc = spsolve(A_all, Rhs)

        # Check convergence
        error = np.amax(np.absolute(Conc - Conc_tmp)) / np.amax(np.absolute(Conc))
        Conc_tmp = Conc.copy()

        # Extract A and B from interleaved solution
        A = Conc[0::2]
        B = Conc[1::2]

    t_iter = time.perf_counter() - t_start_iter

    vars.A = A
    vars.B = B

    return vars, iter_count, error, t_iter


def main() -> None:
    """
    Main driver for Turing pattern simulation.
    """
    print("\n" + "="*70)
    print("TURING PATTERN SIMULATION - DOUBLE DIFFUSION FEM SOLVER")
    print("="*70)

    # Geometry parameters
    geom = GeometryParams(
        x0=0.0,
        y0=0.0,
        lx=5.0,
        ly=5.0
    )

    # Time parameters
    time_params = TimeParams(
        dt=5e-4,
        nt=200
    )

    # Model parameters
    model = ModelParams(
        g_coeff=600.0,  # Kinetic constant
        a_coeff=0.05,   # Growth A
        b_coeff=1.0,    # Growth B
        d_coeff=20.0,   # Diffusivity B
        amp=0.01        # Random noise
    )

    # Solver parameters
    solver = SolverParams(
        tol=0.001,
        iter_max=20
    )

    # Gauss integration points for triangles
    integration = IntegrationParams(
        nip=3,
        gauss=np.array([[1/6, 2/3, 1/6], [1/6, 1/6, 2/3]]).T.copy(),
        weights=np.array([1/6, 1/6, 1/6])
    )

    # Pre-compute shape functions
    print("\nPre-computing shape functions at integration points...")
    integration = precompute_shape_functions(integration)
    print(f"  - Shape functions pre-computed for {integration.nip} integration points")

    # Create mesh
    print("\nGenerating mesh...")
    t_start_mesh = time.perf_counter()
    mesh = make_mesh(geom)
    t_mesh = time.perf_counter() - t_start_mesh
    print(f"Mesh generated in {t_mesh*1000:.2f} ms")
    print(f"  - Number of elements: {mesh.nel}")
    print(f"  - Number of nodes: {mesh.nnod}")
    print(f"  - System DOFs: {mesh.sdof} (2 per node)")

    # Initialize concentrations
    print("\nInitializing concentration fields...")
    vars = initialize_concentrations(mesh, model)
    print(f"  - Species A: mean={np.mean(vars.A):.4f}, std={np.std(vars.A):.4f}")
    print(f"  - Species B: mean={np.mean(vars.B):.4f}, std={np.std(vars.B):.4f}")

    # Setup output writing
    print(f"\nSetting up output files...")
    points = np.hstack((mesh.GCOORD, mesh.GCOORD[:, 0].reshape(-1, 1)*0))
    cells = [("triangle", mesh.EL2NOD)]
    writer = meshio.xdmf.TimeSeriesWriter('transient.xmf')
    writer.__enter__()
    writer.write_points_cells(points, cells)

    # Assemble system matrix (constant during time stepping)
    print("\nAssembling system matrix...")
    t_start_assembly = time.perf_counter()
    I, J, K, Rhs_base = assemble_system_optimized(
        mesh.GCOORD,
        mesh.EL2NOD,
        vars.A,
        vars.B,
        integration.N_all,
        integration.dNds_all,
        integration.weights,
        model.g_coeff,
        model.d_coeff,
        time_params.dt
    )
    A_all = csr_matrix((K.reshape(mesh.nel*2*mesh.nnodel*mesh.nnodel),
                        (I.reshape(mesh.nel*2*mesh.nnodel*mesh.nnodel),
                         J.reshape(mesh.nel*2*mesh.nnodel*mesh.nnodel))),
                       shape=(mesh.sdof, mesh.sdof))
    t_assembly = time.perf_counter() - t_start_assembly
    print(f"System matrix assembled in {t_assembly*1000:.2f} ms")
    print(f"  - Matrix size: {mesh.sdof} x {mesh.sdof}")
    print(f"  - Non-zeros: {A_all.nnz}")

    # Time loop
    print(f"\nStarting time integration ({time_params.nt} steps, dt={time_params.dt})...")
    mode_str = "NUMBA" if (USE_NUMBA and NUMBA_AVAILABLE) else "PYTHON"
    print(f"Execution mode: {mode_str}")
    t_start_timeloop = time.perf_counter()

    total_iterations = 0
    for t in range(time_params.nt):
        # Update base RHS with current concentrations
        _, _, _, Rhs_base = assemble_system_optimized(
            mesh.GCOORD,
            mesh.EL2NOD,
            vars.A,
            vars.B,
            integration.N_all,
            integration.dNds_all,
            integration.weights,
            model.g_coeff,
            model.d_coeff,
            time_params.dt
        )

        # Solve with nonlinear iterations
        vars, iter_count, error, t_iter = solve_timestep(
            mesh, model, time_params, solver, integration,
            vars, A_all, Rhs_base, t
        )

        total_iterations += iter_count

        # Save data
        writer.write_data(t, point_data={"A": vars.A, "B": vars.B})

        # Print progress
        if t == 0 or (t + 1) % 10 == 0:
            print(f"Step {t+1}/{time_params.nt}: "
                  f"iterations={iter_count}, error={error:.2e}, "
                  f"time={t_iter*1000:.2f} ms")

    t_timeloop = time.perf_counter() - t_start_timeloop

    writer.__exit__()

    # Final summary
    print(f"\n{'='*70}")
    print("SIMULATION COMPLETE - PERFORMANCE SUMMARY")
    print(f"{'='*70}")
    print(f"Mesh generation:       {t_mesh*1000:10.2f} ms")
    print(f"Matrix assembly:       {t_assembly*1000:10.2f} ms")
    print(f"Time integration:      {t_timeloop:10.2f} s  ({time_params.nt} steps)")
    print(f"Average per step:      {t_timeloop/time_params.nt*1000:10.2f} ms")
    print(f"Average iterations:    {total_iterations/time_params.nt:.1f}")
    print(f"Total simulation time: {(t_mesh + t_assembly + t_timeloop):10.2f} s")
    print(f"Execution mode:        {mode_str}")
    print(f"{'='*70}")
    print(f"Output files: transient.xmf, transient.h5")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
