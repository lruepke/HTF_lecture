"""
Driver script for 2D mechanical Stokes flow FEM solver.

BLOCKED VECTORIZED version - processes elements in blocks for
better cache utilization and vectorized assembly.
"""

from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
import triangle as tr
from mechanical2d_blocked import (
    Mesh, MaterialParams, BoundaryConditions,
    SolverParams, solve_mechanical_2d
)


@dataclass
class GeometryParams:
    """Container for geometry parameters."""
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    c_inc: list    # List of inclusion centers [(x1, y1), (x2, y2), ...]
    radius: list   # List of radii [r1, r2, ...]
    no_pts: int    # Number of points for inclusion boundary


def make_mesh(geom: GeometryParams, el_ids=(1, 100), el_sizes=(0.1, 0.01)) -> Mesh:
    """
    Create an unstructured triangle mesh with multiple circular inclusions.

    Parameters
    ----------
    geom : GeometryParams
        Geometry parameters including domain size and inclusion properties
    el_ids : tuple
        Material IDs for (matrix, inclusion)
    el_sizes : tuple
        Target element sizes for (matrix, inclusion)

    Returns
    -------
    Mesh
        Mesh object containing node coordinates, connectivity, and phase IDs
    """
    vertices = []
    segments = []
    regions = []
    segment_markers = []

    def _make_box() -> None:
        """Add rectangular domain boundary to mesh."""
        i = len(vertices)
        vertices.extend([
            [geom.x_min, geom.y_min],
            [geom.x_max, geom.y_min],
            [geom.x_max, geom.y_max],
            [geom.x_min, geom.y_max]
        ])
        segments.extend([
            (i+0, i+1), (i+1, i+2), (i+2, i+3), (i+3, i+0)
        ])
        regions.append([
            geom.x_min + 0.01 * (geom.x_max - geom.x_min),
            geom.y_min + 0.01 * (geom.y_max - geom.y_min),
            el_ids[0], el_sizes[0]
        ])
        segment_markers.extend([101, 102, 103, 104])

    def _make_inclusion(center: tuple, radius: float) -> None:
        """Add circular inclusion to mesh.

        Parameters
        ----------
        center : tuple
            Center coordinates (x, y) of the inclusion
        radius : float
            Radius of the inclusion
        """
        theta = np.linspace(0, 2*np.pi, geom.no_pts, endpoint=False)
        xx = np.cos(theta)
        yy = np.sin(theta)

        i = len(vertices)
        vertices.extend(
            np.array([center[0] + radius*xx,
                     center[1] + radius*yy]).T
        )

        Tmp = np.array([np.arange(i, i+geom.no_pts),
                       np.arange(i+1, i+geom.no_pts+1)]).T
        Tmp[-1, 1] = i
        segments.extend(Tmp)
        segment_markers.extend(1001 * np.ones(len(Tmp)))

        regions.append([center[0], center[1], el_ids[1], el_sizes[1]])

    # Build geometry
    _make_box()
    # Add all inclusions
    for center, radius in zip(geom.c_inc, geom.radius):
        _make_inclusion(center, radius)

    # Generate mesh using Triangle
    A = dict(
        vertices=vertices,
        segments=segments,
        segment_markers=segment_markers,
        regions=regions
    )
    B = tr.triangulate(A, 'o2pq33Aa')

    # Extract mesh data
    GCOORD = B.get("vertices")
    EL2NOD = B.get("triangles")
    Phases = B.get("triangle_attributes").flatten().astype(int)
    Node_ids = B.get("vertex_markers").flatten()

    # Add center node to each element (7-node triangle)
    nnod = GCOORD.shape[0]
    GCOORD = np.vstack((GCOORD, np.mean(GCOORD[EL2NOD[:, 0:3]], axis=1)))
    new_nodes = np.arange(nnod, GCOORD.shape[0]).reshape(-1, 1)
    EL2NOD = np.hstack((EL2NOD, new_nodes))
    Node_ids = np.hstack((Node_ids, np.zeros(new_nodes.shape[0]))).astype(int)

    return Mesh(
        GCOORD=GCOORD,
        EL2NOD=EL2NOD,
        Phases=Phases,
        Node_ids=Node_ids,
        nnod=GCOORD.shape[0],
        nel=EL2NOD.shape[0],
        nnodel=EL2NOD.shape[1]
    )


def plot_solution(mesh: Mesh, velocity: np.ndarray, pressure: np.ndarray,
                 filename: str = 'pressure_field.png') -> None:
    """
    Visualize pressure field, mesh, and velocity vectors.

    Parameters
    ----------
    mesh : Mesh
        Mesh object
    velocity : np.ndarray
        Velocity field (2*nnod DOFs: [vx0, vy0, vx1, vy1, ...])
    pressure : np.ndarray
        Pressure field (discontinuous, 3 values per element)
    filename : str
        Output filename for figure
    """
    # Create expanded connectivity for discontinuous pressure visualization
    GC_BIG = np.vstack((
        mesh.GCOORD[mesh.EL2NOD[:, 0:3], 0].ravel(),
        mesh.GCOORD[mesh.EL2NOD[:, 0:3], 1].ravel()
    )).reshape(2, -1).T
    EN_BIG = np.arange(3 * mesh.nel).reshape(-1, 3)

    # Extract velocity components
    Vel_x = velocity[0::2]  # x-components at nodes
    Vel_y = velocity[1::2]  # y-components at nodes

    plt.figure(figsize=(10, 8))
    levels = np.linspace(pressure.min(), pressure.max(), num=100)
    contours = plt.tricontourf(GC_BIG[:, 0], GC_BIG[:, 1], EN_BIG,
                               pressure, levels=levels, cmap='jet')
    plt.triplot(mesh.GCOORD[:, 0], mesh.GCOORD[:, 1], mesh.EL2NOD[:, 0:3],
               color='black', linewidth=0.1, alpha=0.3)

    plt.colorbar(contours, label='Pressure')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.title('Pressure Field and Velocity Vectors')
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Figure saved as {filename}")
    plt.show()


def main() -> None:
    """
    Main driver for 2D mechanical Stokes flow FEM simulation.
    """
    print("\n" + "="*70)
    print("2D STOKES FLOW - FEM SOLVER (BLOCKED VECTORIZED)")
    print("="*70)

    # Geometry parameters with multiple inclusions
    geom = GeometryParams(
        x_min=-1.0,
        x_max=1.0,
        y_min=-1.0,
        y_max=1.0,
        #c_inc=[(0.0, 0.0)],  # Multiple inclusion centers
        #radius=[0.2],                       # Corresponding radii
        c_inc=[(0.0, 0.0), (-0.5, 0.5), (0.5, -0.5), (0.25, 0.25), (0.4,0.05), (0.5, -.15)],  # Multiple inclusion centers
        radius=[0.2, 0.15, 0.18, 0.1, 0.1, 0.1],                       # Corresponding radii
        no_pts=60
    )

    # Material parameters
    material = MaterialParams(
        Mu=np.array([1e0, 2e0]),    # Viscosity [matrix, inclusion]
        Rho=np.array([1.0, 2.0]),   # Density [matrix, inclusion]
        G=np.array([0.0, 0.0])      # Gravity [gx, gy]
    )

    # Solver parameters
    solver = SolverParams(nip=6)

    # Element properties
    el_sizes = (1e-3, 1e-3)  # Target element size [matrix, inclusion]
    el_ids = (0, 1)          # Material IDs [matrix, inclusion]

    # Generate mesh
    print("\nGenerating mesh...")
    mesh = make_mesh(geom=geom, el_ids=el_ids, el_sizes=el_sizes)

    print(f"Mesh generated:")
    print(f"  - Number of nodes: {mesh.nnod}")
    print(f"  - Number of elements: {mesh.nel}")
    print(f"  - Nodes per element: {mesh.nnodel}")

    # Material distribution
    n_matrix = np.sum(mesh.Phases == el_ids[0])
    n_inclusion = np.sum(mesh.Phases == el_ids[1])
    print(f"\nMaterial distribution:")
    print(f"  - Matrix elements (mu={material.Mu[el_ids[0]]}): {n_matrix}")
    print(f"  - Inclusion elements (mu={material.Mu[el_ids[1]]}): {n_inclusion}")

    # Set boundary conditions (pure shear: u = 0.5*x, v = -0.5*y)
    print("\nApplying boundary conditions...")
    Bc_ind = np.where(np.isin(mesh.Node_ids, [101, 102, 103, 104]))[0]
    Bc_val = np.hstack((
        0.5 * mesh.GCOORD[Bc_ind, 0],
        -0.5 * mesh.GCOORD[Bc_ind, 1]
    ))
    Bc_ind = np.hstack((2*Bc_ind, 2*Bc_ind+1))

    bc = BoundaryConditions(Bc_ind=Bc_ind, Bc_val=Bc_val)
    print(f"  - Prescribed DOFs: {len(bc.Bc_ind)}")

    # Solve
    print("\nSolving Stokes flow problem...")
    solution = solve_mechanical_2d(
        mesh=mesh,
        material=material,
        bc=bc,
        solver=solver
    )

    print(f"Solution obtained:")
    print(f"  - Velocity DOFs: {solution.Vel.shape[0]}")
    print(f"  - Pressure DOFs: {solution.Pressure.shape[0]}")
    print(f"  - Max velocity: {np.max(np.abs(solution.Vel)):.6f}")
    print(f"  - Pressure range: [{solution.Pressure.min():.6f}, "
          f"{solution.Pressure.max():.6f}]")

    # Visualize results
    print("\nGenerating visualization...")
    plot_solution(mesh, solution.Vel, solution.Pressure)

    print(f"\n{'='*70}")
    print("SIMULATION COMPLETE")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
