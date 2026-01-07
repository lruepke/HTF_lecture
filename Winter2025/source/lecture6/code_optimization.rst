Code Refactoring and Performance Optimization
==============================================

In this section, we'll explore how to improve our FEM code through refactoring and performance optimization. We'll take the transient heat diffusion solver from :download:`2d_fem_transient_triangle.py <python/2d_fem_transient_triangle.py>` and progressively improve it in two stages:

1. **Refactoring (v2)**: Improving code organization and readability using modern Python features
2. **Optimization (v3)**: Accelerating performance using numba JIT compilation and algorithmic improvements

This demonstrates important software engineering practices: write clear code first, then optimize where needed.

Code Refactoring: From Original to v2
--------------------------------------

Why Refactor?
^^^^^^^^^^^^^

The original implementation works correctly but has several characteristics common in research code:

- Many global or scattered variables
- Function signatures with many positional arguments
- Limited type information
- Implicit data relationships

As projects grow, these patterns make code harder to maintain, debug, and extend. Let's improve this systematically.

Using Dataclasses for Data Organization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Python's :code:`dataclass` decorator (Python 3.7+) provides a clean way to group related data. Instead of passing many individual parameters:

.. code-block:: python
    :caption: Original approach - many parameters

    def solve_fem(GCOORD, EL2NOD, Phases, rho, cp, k1, k2,
                  Ttop, Tbot, dt, nt, nip, gauss, weights, ...):
        # Too many parameters to track!
        pass

We organize data into logical groups:

.. code-block:: python
    :caption: Refactored approach - dataclasses

    from dataclasses import dataclass
    import numpy as np

    @dataclass
    class Mesh:
        """Container for mesh data."""
        GCOORD: np.ndarray  # Node coordinates, shape (nnod, 2)
        EL2NOD: np.ndarray  # Element connectivity, shape (nel, nnodel)
        Phases: np.ndarray  # Phase/material ID per element
        nnod: int = 0       # Number of nodes
        nel: int = 0        # Number of elements
        nnodel: int = 0     # Number of nodes per element

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
    class TimeParams:
        """Container for time stepping parameters."""
        dt: float    # Time step
        nt: int      # Number of time steps

This organization provides several benefits:

- **Self-documenting**: Clear what data belongs together
- **Type hints**: IDEs can provide better autocomplete and error checking
- **Validation**: Can add constraints (e.g., :code:`dt > 0`)
- **Immutability**: Can use :code:`frozen=True` for immutable data
- **Default values**: Easy to specify defaults

Function Signatures with Type Hints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Modern Python supports type hints that make code more readable and enable better tooling:

.. code-block:: python
    :caption: Function with type hints and dataclasses

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
        Solve 2D transient heat diffusion using FEM.

        Parameters
        ----------
        mesh : Mesh
            Mesh object containing coordinates and connectivity
        material : MaterialParams
            Material properties
        bc : BoundaryConditions
            Boundary condition values
        ... [additional parameters]

        Returns
        -------
        Variables
            Updated temperature and conductivity fields
        """
        # Implementation here
        pass

Benefits of this approach:

- **Clear intent**: Immediately see what data is needed
- **Better errors**: Type checkers can catch mistakes before runtime
- **IDE support**: Autocomplete knows the structure of :code:`mesh.GCOORD`, etc.
- **Documentation**: Types serve as inline documentation

Structured Main Function
^^^^^^^^^^^^^^^^^^^^^^^^^

The refactored :code:`main()` function becomes much more readable:

.. code-block:: python
    :caption: Organized main function

    def main() -> None:
        """Main driver for 2D transient heat diffusion FEM simulation."""

        # Define all parameters using dataclasses
        geom = GeometryParams(x0=-1.0, y0=-1.0, lx=2.0, ly=2.0,
                             n_incl=5, radius=0.15)

        time_params = TimeParams(dt=0.015, nt=80)

        material = MaterialParams(rho=1.0, cp=1.0, k1=1.0, k2=0.01)

        bc = BoundaryConditions(Ttop=0.0, Tbot=1.0)

        integration = IntegrationParams(
            nip=3,
            gauss=np.array([[1/6, 2/3, 1/6], [1/6, 1/6, 2/3]]).T,
            weights=np.array([1/6, 1/6, 1/6])
        )

        # Create mesh
        mesh = make_mesh(geom)

        # Initialize variables
        vars = Variables(T=np.zeros(mesh.nnod), Kel=None)

        # Time loop
        for t in range(time_params.nt):
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

The complete refactored code is available here: :download:`2d_fem_transient_triangle_v2.py <python/2d_fem_transient_triangle_v2.py>`.

Performance Optimization: From v2 to v3
----------------------------------------

Performance Profiling
^^^^^^^^^^^^^^^^^^^^^

Before optimizing, we must measure where time is spent. Adding timing instrumentation to v2 reveals:

.. code-block:: python
    :caption: Adding performance timing

    import time

    def solve_2d_temperature_fem(...):
        t_start_total = time.perf_counter()

        # Element assembly
        t_start_assembly = time.perf_counter()
        for iel in range(nel):
            # ... element loop ...
        t_assembly = time.perf_counter() - t_start_assembly

        # ... more timing sections ...

        # Report timings
        print(f"Element assembly: {t_assembly*1000:.2f} ms ({t_assembly/t_total*100:.1f}%)")

Running the profiled v2 code shows:

.. code-block:: text

    Time step 10/80 - Performance breakdown:
    ======================================================================
      Element assembly:    117 ms (75.7%)  ← Primary bottleneck
      Postprocessing:       33 ms (22.4%)  ← Secondary bottleneck
      Sparse matrix:         0.15 ms (0.1%)
      Linear solve:          2.3 ms (1.5%)
      ──────────────────────────────────────
      Total time:          155 ms

**Key finding**: Over 98% of time is spent in Python loops (element assembly and postprocessing). These are prime candidates for optimization.

Optimization Strategy 1: Pre-computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For triangular elements with fixed Gauss integration points, the shape functions :math:`N_i` and their local derivatives :math:`\frac{\partial N_i}{\partial \xi}` are **constant** at each integration point.

Original code computes these 7,635 times per timestep:

.. math::

    \text{Function calls} = n_{\text{elements}} \times n_{\text{integration points}} = 2545 \times 3 = 7635

We can pre-compute them once:

.. code-block:: python
    :caption: Pre-computing shape functions

    def precompute_shape_functions(integration: IntegrationParams) -> IntegrationParams:
        """
        Pre-compute shape functions at all integration points.

        For triangular elements, N and dNds are constant at each
        integration point - no need to recompute every timestep.
        """
        nip = integration.nip
        nnodel = 3

        N_all = np.zeros((nip, nnodel))
        dNds_all = np.zeros((nip, 2, nnodel))

        for ip in range(nip):
            xi = integration.gauss[ip, 0]
            eta = integration.gauss[ip, 1]
            N, dNds = shapes_tri(xi, eta)
            N_all[ip, :] = N
            dNds_all[ip, :, :] = dNds

        integration.N_all = N_all
        integration.dNds_all = dNds_all
        return integration

Usage in element loop:

.. code-block:: python
    :caption: Using pre-computed values

    # OLD: Compute every time
    N, dNds = shapes_tri(xi, eta)

    # NEW: Lookup pre-computed values
    N = integration.N_all[ip, :]
    dNds = integration.dNds_all[ip, :, :]

**Impact**: 5-10% speedup with zero algorithm change—just avoiding redundant computation.

Optimization Strategy 2: Numba JIT Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Numba <https://numba.pydata.org/>`_ is a just-in-time (JIT) compiler that translates Python code to optimized machine code. It excels at accelerating numerical loops.

Setting up Numba with a toggle:

.. code-block:: python
    :caption: Numba configuration at top of file

    # ============================================================================
    # PERFORMANCE CONFIGURATION
    # ============================================================================
    USE_NUMBA = True  # Set to False to disable numba optimization

    from dataclasses import dataclass
    import numpy as np
    # ... other imports ...

    # Try to import numba
    try:
        from numba import njit
        NUMBA_AVAILABLE = True
    except ImportError:
        NUMBA_AVAILABLE = False
        if USE_NUMBA:
            print("Warning: Numba not available, falling back to pure Python")
        USE_NUMBA = False

This provides:

- **Easy toggle**: Students can enable/disable with one line
- **Graceful fallback**: Works without numba installation
- **Educational value**: Can compare implementations side-by-side

Creating Numba-Optimized Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We create numba-compiled versions of our computational kernels:

.. code-block:: python
    :caption: Numba-optimized element assembly

    @njit(cache=True)
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
        """Numba-optimized element assembly loop."""
        nel = EL2NOD.shape[0]
        nnodel = EL2NOD.shape[1]
        nnod = GCOORD.shape[0]
        nip = N_all.shape[0]

        # Storage arrays
        Rhs_all = np.zeros(nnod)
        I = np.zeros((nel, nnodel * nnodel))
        J = np.zeros((nel, nnodel * nnodel))
        K = np.zeros((nel, nnodel * nnodel))

        for iel in range(nel):
            ECOORD = GCOORD[EL2NOD[iel, :], :]
            Ael = np.zeros((nnodel, nnodel))
            Rhs_el = np.zeros(nnodel)

            for ip in range(nip):
                N = N_all[ip, :]
                dNds = dNds_all[ip, :, :]

                # Jacobian
                Jac = dNds @ ECOORD
                invJ = np.linalg.inv(Jac)
                detJ = np.linalg.det(Jac)
                dNdx = invJ @ dNds

                # Element matrices
                Me = np.outer(N, N)
                Ke_local = dNdx.T @ dNdx

                Ael += (rho * cp * Me + dt * Kel[iel] * Ke_local) * detJ * weights[ip]

                T_el = T[EL2NOD[iel, :]]
                Rhs_el += rho * cp * (Me @ T_el) * detJ * weights[ip]

            # Store sparse matrix entries
            for i in range(nnodel):
                for j in range(nnodel):
                    idx = i * nnodel + j
                    I[iel, idx] = EL2NOD[iel, i]
                    J[iel, idx] = EL2NOD[iel, j]
                    K[iel, idx] = Ael[i, j]

            for i in range(nnodel):
                Rhs_all[EL2NOD[iel, i]] += Rhs_el[i]

        return I, J, K, Rhs_all

The :code:`@njit(cache=True)` decorator:

- Compiles the function to machine code on first call (slow first run)
- Caches compiled version for subsequent runs (fast afterward)
- Enforces type consistency for optimization

Conditional Dispatch
^^^^^^^^^^^^^^^^^^^^

We keep both Python and Numba code paths for comparison:

.. code-block:: python
    :caption: Dual implementation paths

    def solve_2d_temperature_fem(...):
        # ... setup code ...

        if USE_NUMBA and NUMBA_AVAILABLE:
            # ========== NUMBA PATH ==========
            I, J, K, Rhs_all = assemble_system_numba(
                mesh.GCOORD, mesh.EL2NOD, Kel, T,
                integration.N_all, integration.dNds_all,
                integration.weights, material.rho, material.cp,
                time_params.dt
            )
        else:
            # ========== PYTHON PATH ==========
            # Original implementation (unchanged)
            for iel in range(nel):
                # ... element loop ...

        # ... rest of solver ...

This design:

- Preserves original code for teaching/debugging
- Allows easy performance comparison
- Students can see both implementations

Optimization Strategy 3: On-Demand Computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Profiling the numba-optimized code reveals a surprising result:

.. code-block:: text

    Time step 10/80 - Performance breakdown (Mode: NUMBA):
    ======================================================================
      Element assembly:     8.5 ms (35.1%)
      Postprocessing:      12.8 ms (52.3%)  ← Now the bottleneck!
      Linear solve:         2.5 ms (10.4%)
      ──────────────────────────────────────
      Total time:          24.4 ms

**Problem**: Postprocessing (computing heat flux for visualization) takes more time than assembly! But we only write output every 10 timesteps.

**Solution**: Compute flux only when actually writing output:

.. code-block:: python
    :caption: On-demand postprocessing

    @dataclass
    class TimeParams:
        dt: float
        nt: int
        output_freq: int  # Write output every N timesteps

    def solve_2d_temperature_fem(...):
        # ... solve system ...

        # Check if we should write output
        write_output = (t % time_params.output_freq == 0) or (t == time_params.nt - 1)

        if write_output:
            # Only compute flux when needed
            Q_x, Q_y = compute_heat_flux_numba(...)
            writer.write_data(t, point_data={"T": T},
                            cell_data={"U": [U], "K": [Kel]})
        else:
            # Skip postprocessing entirely
            pass

This algorithmic optimization provides massive speedup:

.. code-block:: text

    Non-output steps (90% of timesteps):
      Element assembly:     7.5 ms (76.5%)
      Linear solve:         2.1 ms (21.4%)
      Postprocessing:       SKIPPED
      ──────────────────────────────────────
      Total time:           9.8 ms

**Lesson**: Sometimes the best optimization is *not computing* something at all!

Performance Results
-------------------

Summary of Optimization Impact
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Performance Evolution
   :widths: 30 20 20 30
   :header-rows: 1

   * - Version
     - Time/step
     - Speedup
     - Key Features
   * - Original Python
     - 159 ms
     - 1.0×
     - Baseline implementation
   * - v2 (Refactored)
     - 155 ms
     - 1.03×
     - Better code organization
   * - v3 (Pre-compute)
     - 145 ms
     - 1.10×
     - Shape function caching
   * - v3 (+ Numba)
     - 24 ms
     - 6.6×
     - JIT compilation
   * - v3 (+ On-demand)
     - 11 ms
     - 14.5×
     - Skip unnecessary work

Detailed Timing Breakdown
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: _static/figures/performance_comparison.png
   :width: 80%
   :align: center

   Time distribution across different implementations. Note how bottlenecks shift: assembly dominates initially, then postprocessing becomes critical after numba optimization, finally resolved by on-demand computation.

The cumulative effect is dramatic:

.. math::

    \text{Total simulation time (80 steps)}:\ \ 12.7\,\text{s} \rightarrow 0.9\,\text{s}

Educational Value
-----------------

Code Quality vs Performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This progression teaches important lessons:

1. **Refactor first, optimize second**

   - v2 improves maintainability with no performance cost
   - Clean code is easier to optimize later
   - Type hints and dataclasses help catch bugs

2. **Measure before optimizing**

   - Timing revealed postprocessing bottleneck
   - Intuition about bottlenecks is often wrong
   - Profile-guided optimization is essential

3. **Multiple optimization strategies**

   - Algorithmic (pre-computation, on-demand)
   - Compilation (numba)
   - Each has different impact and complexity

4. **Trade-offs matter**

   - On-demand: less frequent output vs faster simulation
   - Numba: compilation time vs runtime speed
   - Code complexity vs performance gain

Choosing When to Optimize
^^^^^^^^^^^^^^^^^^^^^^^^^^

Not all code needs optimization. Consider:

.. list-table:: When to Optimize
   :widths: 40 60
   :header-rows: 1

   * - Situation
     - Recommendation
   * - Exploratory research
     - Prioritize code clarity (v2)
   * - Production runs
     - Full optimization (v3)
   * - Teaching/debugging
     - Keep both paths, toggle as needed
   * - Quick prototyping
     - Start with v2, optimize if slow

The complete optimized code is available here: :download:`2d_fem_transient_triangle_v3.py <python/2d_fem_transient_triangle_v3.py>`.

Exercises
---------

1. **Profiling practice**

   Add timing to your own FEM code. What fraction of time is spent in each section? Are the bottlenecks where you expected?

2. **Numba exploration**

   Modify :code:`USE_NUMBA` to compare Python vs numba performance. Run multiple times to see the impact of JIT compilation caching.

3. **Output frequency experiment**

   Change :code:`output_freq` in v3. How does simulation time scale with output frequency? Plot the relationship.

4. **Code refactoring**

   Take your steady-state FEM code and refactor it using dataclasses. Does it improve readability?

5. **Further optimization**

   The element assembly loop could be parallelized using :code:`numba.prange`. What speedup would you expect on a 4-core CPU? (Hint: Amdahl's law)

Additional Resources
--------------------

- `Python dataclasses documentation <https://docs.python.org/3/library/dataclasses.html>`_
- `Numba user guide <https://numba.readthedocs.io/>`_
- `Performance profiling with cProfile <https://docs.python.org/3/library/profile.html>`_
- :download:`Performance analysis summary <python/PERFORMANCE_SUMMARY.md>`
- :download:`Optimization strategy guide <python/OPTIMIZATION_STRATEGY.md>`
