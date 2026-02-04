"""
Driver script for viscous single-layer folding simulation.

A competent layer (high viscosity) embedded in a weak matrix develops
folds under pure shear compression. Uses the blocked Stokes FEM solver
with Lagrangian mesh advection.

Physics:
    - Pure shear BCs: vx = -eps_bg * x, vy = +eps_bg * y
    - Lagrangian advection: mesh nodes move with the velocity field
    - Dominant wavelength: lambda_d ~ 2*pi*h*(R/6)^(1/3)
"""

import numpy as np
import triangle as tr
import h5py
from folding_solver import (
    Mesh, MaterialParams, BoundaryConditions, SolverParams,
    FoldingParams, TimeParams,
    solve_mechanical_2d, advect_mesh, update_boundary_conditions,
)


def make_folding_mesh(fp: FoldingParams) -> Mesh:
    """
    Create an unstructured triangle mesh with a folding layer.

    The domain is rectangular with a horizontal layer defined by two
    sinusoidal interfaces. Three regions are created:
      - matrix below the layer  (phase 0)
      - the competent layer     (phase 1)
      - matrix above the layer  (phase 0)

    Parameters
    ----------
    fp : FoldingParams
        Folding parameters (domain size, layer geometry, element sizes)

    Returns
    -------
    Mesh
        Mesh object with 7-node triangular elements
    """
    x_min = -fp.Lx / 2
    x_max =  fp.Lx / 2
    y_min = -fp.Ly / 2
    y_max =  fp.Ly / 2

    vertices = []
    segments = []
    segment_markers = []
    regions = []

    # ---- Domain boundary (4 corners) ----
    i0 = len(vertices)
    vertices.extend([
        [x_min, y_min],   # 0 - bottom-left
        [x_max, y_min],   # 1 - bottom-right
        [x_max, y_max],   # 2 - top-right
        [x_min, y_max],   # 3 - top-left
    ])
    # Bottom edge (marker 101), right (102), top (103), left (104)
    segments.extend([
        (i0 + 0, i0 + 1),
        (i0 + 1, i0 + 2),
        (i0 + 2, i0 + 3),
        (i0 + 3, i0 + 0),
    ])
    segment_markers.extend([101, 102, 103, 104])

    # ---- Discretize bottom interface ----
    n_pts = fp.n_interface_pts
    x_pts = np.linspace(x_min, x_max, n_pts)
    y_bot = -fp.h / 2 + fp.A * np.cos(2 * np.pi * x_pts / fp.wavelength)
    y_top =  fp.h / 2 + fp.A * np.cos(2 * np.pi * x_pts / fp.wavelength)

    # Add bottom interface vertices
    i_bot_start = len(vertices)
    for k in range(n_pts):
        vertices.append([x_pts[k], y_bot[k]])
    # Segments along the bottom interface (marker 201)
    for k in range(n_pts - 1):
        segments.append((i_bot_start + k, i_bot_start + k + 1))
        segment_markers.append(201)

    # ---- Add top interface vertices ----
    i_top_start = len(vertices)
    for k in range(n_pts):
        vertices.append([x_pts[k], y_top[k]])
    # Segments along the top interface (marker 202)
    for k in range(n_pts - 1):
        segments.append((i_top_start + k, i_top_start + k + 1))
        segment_markers.append(202)

    # ---- Connect interfaces to domain boundary at left and right edges ----
    # Left side connections
    segments.append((i0 + 3, i_top_start))
    segment_markers.append(104)
    segments.append((i_top_start, i_bot_start))
    segment_markers.append(104)
    segments.append((i_bot_start, i0 + 0))
    segment_markers.append(104)

    # Right side connections
    segments.append((i0 + 1, i_bot_start + n_pts - 1))
    segment_markers.append(102)
    segments.append((i_bot_start + n_pts - 1, i_top_start + n_pts - 1))
    segment_markers.append(102)
    segments.append((i_top_start + n_pts - 1, i0 + 2))
    segment_markers.append(102)

    # Remove original left (index 3) and right (index 1) boundary segments
    # since they are now split by the interface connections above
    new_segments = []
    new_markers = []
    for idx, (seg, mk) in enumerate(zip(segments, segment_markers)):
        if idx == 1 or idx == 3:
            continue
        new_segments.append(seg)
        new_markers.append(mk)
    segments = new_segments
    segment_markers = new_markers

    # ---- Region points ----
    y_region_below = (y_min + y_bot.min()) / 2
    regions.append([0.0, y_region_below, 0, fp.el_size_matrix])
    regions.append([0.0, 0.0, 1, fp.el_size_layer])
    y_region_above = (y_max + y_top.max()) / 2
    regions.append([0.0, y_region_above, 0, fp.el_size_matrix])

    # ---- Triangulate ----
    A = dict(
        vertices=vertices,
        segments=segments,
        segment_markers=segment_markers,
        regions=regions,
    )
    B = tr.triangulate(A, 'o2pq33Aa')

    GCOORD = B.get("vertices")
    EL2NOD = B.get("triangles")
    Phases = B.get("triangle_attributes").flatten().astype(int)
    Node_ids = B.get("vertex_markers").flatten()

    # ---- Fix boundary markers ----
    # Interface endpoint nodes on the domain boundary may receive interface
    # markers (201/202) from triangle instead of boundary markers. Fix by
    # coordinate: any node on the domain edge must have a boundary marker.
    tol = 1e-10 * max(fp.Lx, fp.Ly)
    on_bottom = np.abs(GCOORD[:, 1] - y_min) < tol
    on_right  = np.abs(GCOORD[:, 0] - x_max) < tol
    on_top    = np.abs(GCOORD[:, 1] - y_max) < tol
    on_left   = np.abs(GCOORD[:, 0] - x_min) < tol

    not_boundary = ~np.isin(Node_ids, [101, 102, 103, 104])
    Node_ids[on_bottom & not_boundary] = 101
    Node_ids[on_right  & not_boundary] = 102
    Node_ids[on_top    & not_boundary] = 103
    Node_ids[on_left   & not_boundary] = 104

    # Add 7th (center) node to each element
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
        nnodel=EL2NOD.shape[1],
    )


# ============================================================================
# XDMF/HDF5 OUTPUT (manual, supports deforming mesh)
# ============================================================================

def write_output_step(h5file, step_info, mesh, material, velocity, pressure,
                      time_val):
    """
    Write one timestep to the HDF5 file.

    Uses an expanded mesh (3 unique corner nodes per element) so the
    discontinuous P-1 pressure is written as smooth point data.

    Parameters
    ----------
    h5file : h5py.File
        Open HDF5 file
    step_info : list
        Accumulates step metadata (modified in-place)
    mesh : Mesh
        Current (possibly deformed) mesh
    material : MaterialParams
        Material properties (for per-element viscosity)
    velocity : np.ndarray
        Velocity field (2*nnod,)
    pressure : np.ndarray
        Discontinuous pressure (3*nel,)
    time_val : float
        Simulation time
    """
    nel = mesh.nel
    corners = mesh.EL2NOD[:, 0:3]  # (nel, 3)

    # Expanded points: each element gets its own 3 corner nodes
    points = np.column_stack([
        mesh.GCOORD[corners, 0].ravel(),
        mesh.GCOORD[corners, 1].ravel(),
        np.zeros(3 * nel),
    ])

    # Sequential connectivity (constant, but stored per step for generality)
    cells = np.arange(3 * nel, dtype=np.int32).reshape(-1, 3)

    # Velocity at expanded corner nodes
    Vx = velocity[0::2]
    Vy = velocity[1::2]
    vel_3d = np.column_stack([
        Vx[corners.ravel()],
        Vy[corners.ravel()],
        np.zeros(3 * nel),
    ])

    # Write datasets
    step_name = f"step_{len(step_info):04d}"
    grp = h5file.create_group(step_name)
    grp.create_dataset("points", data=points)
    grp.create_dataset("cells", data=cells)
    grp.create_dataset("velocity", data=vel_3d)
    grp.create_dataset("pressure", data=pressure.astype(np.float64))
    grp.create_dataset("phase", data=mesh.Phases.astype(np.float64))
    grp.create_dataset("viscosity", data=material.Mu[mesh.Phases])

    step_info.append({
        "name": step_name,
        "time": time_val,
        "n_points": 3 * nel,
        "n_cells": nel,
    })
    print(f"  Written t = {time_val:.4f}")


def finalize_output(xmf_filename, h5_filename, step_info):
    """
    Write the XDMF XML descriptor for the timeseries.

    Each timestep has its own Geometry and Topology so ParaView
    shows the deforming mesh.

    Parameters
    ----------
    xmf_filename : str
        Output .xmf path
    h5_filename : str
        Companion .h5 path (same basename)
    step_info : list of dict
        Step metadata from write_output_step
    """
    h5_base = h5_filename.split("/")[-1]

    lines = [
        '<?xml version="1.0"?>',
        '<Xdmf Version="3.0">',
        '  <Domain>',
        '    <Grid Name="TimeSeries" GridType="Collection"'
        ' CollectionType="Temporal">',
    ]

    for info in step_info:
        n = info["name"]
        np_ = info["n_points"]
        nc = info["n_cells"]
        t = info["time"]
        lines.extend([
            f'      <Grid Name="{n}" GridType="Uniform">',
            f'        <Time Value="{t}"/>',
            f'        <Geometry GeometryType="XYZ">',
            f'          <DataItem Dimensions="{np_} 3"'
            f' NumberType="Float" Precision="8" Format="HDF">',
            f'            {h5_base}:/{n}/points',
            f'          </DataItem>',
            f'        </Geometry>',
            f'        <Topology TopologyType="Triangle"'
            f' NumberOfElements="{nc}">',
            f'          <DataItem Dimensions="{nc} 3"'
            f' NumberType="Int" Precision="4" Format="HDF">',
            f'            {h5_base}:/{n}/cells',
            f'          </DataItem>',
            f'        </Topology>',
            f'        <Attribute Name="Velocity" AttributeType="Vector"'
            f' Center="Node">',
            f'          <DataItem Dimensions="{np_} 3"'
            f' NumberType="Float" Precision="8" Format="HDF">',
            f'            {h5_base}:/{n}/velocity',
            f'          </DataItem>',
            f'        </Attribute>',
            f'        <Attribute Name="Pressure" AttributeType="Scalar"'
            f' Center="Node">',
            f'          <DataItem Dimensions="{np_}"'
            f' NumberType="Float" Precision="8" Format="HDF">',
            f'            {h5_base}:/{n}/pressure',
            f'          </DataItem>',
            f'        </Attribute>',
            f'        <Attribute Name="Phase" AttributeType="Scalar"'
            f' Center="Cell">',
            f'          <DataItem Dimensions="{nc}"'
            f' NumberType="Float" Precision="8" Format="HDF">',
            f'            {h5_base}:/{n}/phase',
            f'          </DataItem>',
            f'        </Attribute>',
            f'        <Attribute Name="Viscosity" AttributeType="Scalar"'
            f' Center="Cell">',
            f'          <DataItem Dimensions="{nc}"'
            f' NumberType="Float" Precision="8" Format="HDF">',
            f'            {h5_base}:/{n}/viscosity',
            f'          </DataItem>',
            f'        </Attribute>',
            f'      </Grid>',
        ])

    lines.extend([
        '    </Grid>',
        '  </Domain>',
        '</Xdmf>',
    ])

    with open(xmf_filename, "w") as f:
        f.write("\n".join(lines))


# ============================================================================
# MAIN DRIVER
# ============================================================================

def main() -> None:
    """
    Main driver for viscous single-layer folding simulation.
    """
    print("\n" + "=" * 70)
    print("VISCOUS SINGLE-LAYER FOLDING SIMULATION")
    print("=" * 70)

    # ---- Parameters ----
    fp = FoldingParams(
        Lx=3.0,
        Ly=3.0,
        h=0.2,
        A=0.02,
        wavelength=3.0,       # = Lx, one full wave
        mu_matrix=1.0,
        mu_layer=100.0,
        eps_bg=0.5,
        el_size_layer=0.02,
        el_size_matrix=0.05,
        n_interface_pts=100,
    )

    tp = TimeParams(
        dt=0.005,
        nt=100,
        output_freq=5,
    )

    R = fp.mu_layer / fp.mu_matrix
    lambda_d = 2 * np.pi * fp.h * (R / 6) ** (1.0 / 3)

    print(f"\nFolding parameters:")
    print(f"  Domain:            {fp.Lx} x {fp.Ly}")
    print(f"  Layer thickness:   {fp.h}")
    print(f"  Perturbation:      A = {fp.A} (A/h = {fp.A/fp.h:.1%})")
    print(f"  Wavelength:        {fp.wavelength}")
    print(f"  Viscosity contrast: R = {R:.0f}")
    print(f"  Dominant wavelength (theory): {lambda_d:.3f}")
    print(f"  Background strain rate: {fp.eps_bg}")
    print(f"\nTime stepping:")
    print(f"  dt = {tp.dt}, nt = {tp.nt}, output every {tp.output_freq} steps")
    print(f"  Total strain: {fp.eps_bg * tp.dt * tp.nt:.3f}")

    # ---- Material ----
    material = MaterialParams(
        Mu=np.array([fp.mu_matrix, fp.mu_layer]),
        Rho=np.array([0.0, 0.0]),   # No gravity-driven flow
        G=np.array([0.0, 0.0]),
    )

    solver = SolverParams(nip=6)

    # ---- Generate mesh ----
    print("\nGenerating mesh...")
    mesh = make_folding_mesh(fp)
    print(f"  Nodes: {mesh.nnod}, Elements: {mesh.nel}")
    n_layer = np.sum(mesh.Phases == 1)
    n_matrix = np.sum(mesh.Phases == 0)
    print(f"  Matrix elements: {n_matrix}, Layer elements: {n_layer}")

    # ---- Initial boundary conditions ----
    bc = update_boundary_conditions(mesh, fp.eps_bg)
    print(f"  Boundary DOFs: {len(bc.Bc_ind)}")

    # ---- Time loop ----
    print(f"\n{'='*70}")
    print("STARTING TIME LOOP")
    print(f"{'='*70}")

    xmf_fname = "folding.xmf"
    h5_fname = "folding.h5"

    h5file = h5py.File(h5_fname, "w")
    step_info = []

    for step in range(tp.nt):
        t_current = step * tp.dt
        print(f"\n--- Step {step}/{tp.nt}, t = {t_current:.4f} ---")

        # 1. Solve Stokes
        solution = solve_mechanical_2d(mesh, material, bc, solver)

        # 2. Output (before advection, so mesh matches the solved velocity)
        if step % tp.output_freq == 0:
            write_output_step(h5file, step_info, mesh, material,
                              solution.Vel, solution.Pressure, t_current)

        # 3. Advect mesh
        advect_mesh(mesh, solution.Vel, tp.dt)

        # 4. Update boundary conditions on deformed mesh
        bc = update_boundary_conditions(mesh, fp.eps_bg)

    # ---- Final output ----
    t_final = tp.nt * tp.dt
    print(f"\n--- Final step, t = {t_final:.4f} ---")
    solution = solve_mechanical_2d(mesh, material, bc, solver)
    write_output_step(h5file, step_info, mesh, material,
                      solution.Vel, solution.Pressure, t_final)

    h5file.close()
    finalize_output(xmf_fname, h5_fname, step_info)

    print(f"\n{'='*70}")
    print("SIMULATION COMPLETE")
    print(f"  Total time steps: {tp.nt}")
    print(f"  Total strain: {fp.eps_bg * t_final:.3f}")
    print(f"  Output: {xmf_fname}, {h5_fname}")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
