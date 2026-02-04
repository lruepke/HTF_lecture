.. _viscous-folding:

Viscous single-layer folding
============================

In this section we apply the Stokes solver to a classic problem in structural geology: the development of folds in a competent layer embedded in a weak matrix under compression. This is a time-dependent problem that requires Lagrangian mesh advection on top of the Stokes solver we have already developed.


Physical setup
--------------

A thin horizontal layer of high viscosity :math:`\eta_\text{layer}` is embedded in an infinite matrix of low viscosity :math:`\eta_\text{matrix}`. The system is subjected to horizontal compression (pure shear). If the layer has a small initial perturbation, the perturbation grows exponentially with time because the competent layer resists shortening and instead buckles.

The boundary conditions are pure shear:

.. math::

    v_x = -\dot{\epsilon}_\text{bg} \, x, \qquad v_y = +\dot{\epsilon}_\text{bg} \, y

where :math:`\dot{\epsilon}_\text{bg}` is the background strain rate. Horizontal compression (:math:`v_x`) is balanced by vertical extension (:math:`v_y`) to maintain incompressibility.


Dominant wavelength
^^^^^^^^^^^^^^^^^^^

Linear stability analysis of this problem :cite:`Biot1961,Ramberg1963` predicts that one wavelength grows faster than all others. This *dominant wavelength* depends on the layer thickness :math:`h` and the viscosity contrast :math:`R = \eta_\text{layer} / \eta_\text{matrix}`:

.. math::

    \lambda_d \approx 2\pi h \left(\frac{R}{6}\right)^{1/3}

For example, with :math:`h = 0.2` and :math:`R = 100`, the dominant wavelength is :math:`\lambda_d \approx 3.2`. In our simulation the domain width is chosen to approximately match one dominant wavelength, so a single fold develops.


Lagrangian advection
--------------------

The Stokes equation gives us an *instantaneous* velocity field for the current geometry. To evolve the fold in time, we use **Lagrangian advection**: at each time step, we move the mesh nodes with the computed velocity field.

.. math::

    \mathbf{x}^{n+1} = \mathbf{x}^n + \mathbf{v}^n \, \Delta t

This is a simple forward Euler integration. The mesh deforms with the material, so the layer interface is always resolved by element boundaries. No remeshing is performed; elements deform along with the material. This works well as long as elements do not become too distorted (which limits the total strain we can simulate).


Time-stepping algorithm
^^^^^^^^^^^^^^^^^^^^^^^

The time loop follows these steps:

1. **Solve Stokes**: compute the velocity and pressure fields on the current mesh.
2. **Output**: write the current state (mesh, velocity, pressure) to file.
3. **Advect mesh**: move all nodes by :math:`\mathbf{v} \Delta t`.
4. **Update center nodes**: recompute the 7th (bubble) node of each element as the centroid of the three corner nodes.
5. **Update boundary conditions**: recompute the pure shear velocities at the new boundary node positions.
6. Repeat.

In code, the advection function is:

.. code-block:: python

    def advect_mesh(mesh, velocity, dt):
        Vx = velocity[0::2]
        Vy = velocity[1::2]

        # Number of real nodes (corners + edge midpoints)
        n_real = mesh.nnod - mesh.nel

        # Move real nodes
        mesh.GCOORD[:n_real, 0] += Vx[:n_real] * dt
        mesh.GCOORD[:n_real, 1] += Vy[:n_real] * dt

        # Recompute center nodes as mean of corner nodes
        mesh.GCOORD[n_real:] = np.mean(
            mesh.GCOORD[mesh.EL2NOD[:, 0:3]], axis=1
        )

Center (bubble) nodes do not carry independent velocity DOFs that should be advected freely; they are constrained to remain at the element centroid, so we recompute them from the (already advected) corner nodes.


Mesh generation
---------------

The mesh must conform to the layer interfaces so that the phase boundary coincides with element edges. We define the top and bottom interfaces of the layer as sinusoidal curves:

.. math::

    y_\text{top}(x) &= +\frac{h}{2} + A \cos\!\left(\frac{2\pi x}{\lambda}\right)

    y_\text{bot}(x) &= -\frac{h}{2} + A \cos\!\left(\frac{2\pi x}{\lambda}\right)

where :math:`A` is the initial perturbation amplitude (typically 10% of the layer thickness) and :math:`\lambda` is the perturbation wavelength. Each interface is discretized into ``n_interface_pts`` points that are added as constrained segments in the Triangle mesh generator.

This divides the domain into three regions:

- **Matrix below** the layer (phase 0)
- **Competent layer** (phase 1)
- **Matrix above** the layer (phase 0)

Each region has its own target element size. Finer elements are used inside and near the layer (:math:`\sim 0.02`) to resolve the folding, while coarser elements are used in the far-field matrix (:math:`\sim 0.05`).

.. code-block:: python

    # Discretize interfaces
    x_pts = np.linspace(x_min, x_max, n_interface_pts)
    y_bot = -h/2 + A * np.cos(2*np.pi*x_pts/wavelength)
    y_top =  h/2 + A * np.cos(2*np.pi*x_pts/wavelength)

    # Add as constrained segments to triangle input
    # Region points placed in each of the 3 zones
    # Triangle generates conforming mesh with 'o2pq33Aa' options


Boundary marker fix
^^^^^^^^^^^^^^^^^^^

A subtlety arises where the layer interfaces meet the lateral domain boundaries. The interface endpoint nodes lie on the domain edge but may receive interface markers (201/202) from the Triangle library instead of boundary markers (101--104). If this happens, those nodes are missing from the boundary condition set and the lateral boundary develops kinks.

The fix is to post-process the markers after mesh generation: any node lying on the domain boundary (identified by coordinate) is reassigned to the correct boundary marker.

.. code-block:: python

    tol = 1e-10 * max(Lx, Ly)
    on_left  = np.abs(GCOORD[:, 0] - x_min) < tol
    on_right = np.abs(GCOORD[:, 0] - x_max) < tol
    # ... same for top/bottom

    not_boundary = ~np.isin(Node_ids, [101, 102, 103, 104])
    Node_ids[on_left  & not_boundary] = 104
    Node_ids[on_right & not_boundary] = 102
    # etc.


Updating boundary conditions
-----------------------------

After each advection step, the boundary nodes have moved to new positions. The pure shear boundary conditions must be recomputed at the new coordinates:

.. code-block:: python

    def update_boundary_conditions(mesh, eps_bg):
        bc_nodes = np.where(
            np.isin(mesh.Node_ids, [101, 102, 103, 104])
        )[0]

        Bc_val_x = -eps_bg * mesh.GCOORD[bc_nodes, 0]
        Bc_val_y =  eps_bg * mesh.GCOORD[bc_nodes, 1]

        Bc_val = np.hstack((Bc_val_x, Bc_val_y))
        Bc_ind = np.hstack((2*bc_nodes, 2*bc_nodes+1))
        return BoundaryConditions(Bc_ind=Bc_ind, Bc_val=Bc_val)

Note that the boundary node identification uses markers (not coordinates), which is important because the coordinates change every time step.


Output format
-------------

The simulation writes a single XDMF/HDF5 file pair (``folding.xmf`` + ``folding.h5``) that can be opened directly in ParaView. Since the mesh deforms over time, the XDMF file contains a separate ``<Geometry>`` and ``<Topology>`` block for each time step, so ParaView shows the evolving mesh.

To visualize the discontinuous pressure field smoothly, an *expanded connectivity* is used: each element gets its own copy of 3 corner nodes (no sharing between elements). The 3 pressure DOFs per element then map directly to nodal point data on this expanded mesh, producing smooth contours within each element.

The output includes:

- **Velocity** (point data, vector): the velocity field at each node.
- **Pressure** (point data, scalar): the discontinuous P-1 pressure on the expanded mesh.
- **Phase** (cell data, scalar): 0 for matrix, 1 for layer.
- **Viscosity** (cell data, scalar): the element viscosity.


Default parameters
-------------------

The default simulation parameters are:

.. list-table::
   :widths: 40 20 40
   :header-rows: 1

   * - Parameter
     - Value
     - Description
   * - :math:`L_x \times L_y`
     - :math:`3 \times 3`
     - Domain size
   * - :math:`h`
     - 0.2
     - Layer thickness
   * - :math:`A`
     - 0.02
     - Perturbation amplitude (10% of :math:`h`)
   * - :math:`\lambda`
     - 3.0
     - Perturbation wavelength (= :math:`L_x`)
   * - :math:`\eta_\text{matrix}`
     - 1
     - Matrix viscosity
   * - :math:`\eta_\text{layer}`
     - 100
     - Layer viscosity (:math:`R = 100`)
   * - :math:`\dot{\epsilon}_\text{bg}`
     - 0.5
     - Background strain rate
   * - :math:`\Delta t`
     - 0.005
     - Time step
   * - :math:`n_t`
     - 100
     - Number of time steps
   * - output frequency
     - 5
     - Write every 5th step


Running the simulation
-----------------------

Run the simulation from the ``python/`` directory:

.. code-block:: bash

    python folding_driver.py

This generates the output files ``folding.xmf`` and ``folding.h5``. Open ``folding.xmf`` in ParaView to visualize the results. Use the time controls to step through the simulation and observe the fold growing.

You can download the code here: :download:`folding_solver.py <python/folding_solver.py>` and :download:`folding_driver.py <python/folding_driver.py>`. The solver also requires :download:`ip_triangle.py <python/ip_triangle.py>` and :download:`shp_deriv_triangle.py <python/shp_deriv_triangle.py>`.


Exercises
---------

1. Run the default simulation and verify that a single fold develops. Open the output in ParaView and color by Phase, Pressure, and Velocity magnitude.

2. **Viscosity contrast**: Change the viscosity contrast :math:`R` (e.g., :math:`R = 10, 50, 200`) and observe how the fold shape and growth rate change. Compare the fold wavelength to the theoretical dominant wavelength :math:`\lambda_d`.

3. **Perturbation wavelength**: Set the perturbation wavelength to half the domain width (:math:`\lambda = L_x / 2`) so that two folds fit in the domain. Does the number of folds match the initial perturbation, or does the dominant wavelength take over?

4. **Time step and total strain**: Increase the number of time steps or the time step size to reach larger total strains. At what point do elements become too distorted for the solver to work? What is the maximum fold amplitude you can reach without remeshing?

5. **Amplitude growth**: Track the fold amplitude (maximum :math:`y`-coordinate of the top interface) as a function of time. In the linear regime, the amplitude should grow exponentially. Plot :math:`\ln(A/A_0)` vs. strain and measure the growth rate. Compare with the analytical amplification factor from thin-plate theory.
