"""
Validation and Stability Checks for FEM Simulations

This module provides utilities for checking numerical stability, validating
input parameters, and monitoring solver convergence. These checks help prevent
silent failures and ensure physically meaningful results.

Author: HTF Lecture Series
Date: 2025
"""

import numpy as np
import warnings


def check_cfl_condition(dt, dx, diffusivity, method='implicit', ndim=2):
    """
    Check CFL (Courant-Friedrichs-Lewy) stability condition for diffusion.

    The CFL condition ensures numerical stability for time-dependent PDEs.
    For diffusion equations, it relates the time step to spatial resolution
    and diffusion coefficient.

    Parameters
    ----------
    dt : float
        Time step [s]
    dx : float
        Characteristic mesh spacing (minimum element size) [m]
    diffusivity : float
        Diffusion coefficient [m²/s]
        For thermal problems: diffusivity = k/(rho*cp)
        For species diffusion: diffusivity = D
    method : str, optional
        Time integration method. Default is 'implicit'.
        Options:
            'explicit' - Forward Euler, explicit FDM/FEM
            'implicit' - Backward Euler, Crank-Nicolson
    ndim : int, optional
        Number of spatial dimensions (1, 2, or 3). Default is 2.

    Returns
    -------
    is_stable : bool
        True if stability condition is satisfied
    cfl_number : float
        Computed CFL number (dimensionless)
    cfl_max : float
        Maximum allowed CFL number for stability
        Returns np.inf for unconditionally stable implicit methods

    Notes
    -----
    **CFL Condition for Diffusion:**

    For explicit methods:
        CFL = D * dt / dx² < 1/(2*ndim)

    where ndim is the number of spatial dimensions.

    Stability limits:
        - 1D: dt < dx² / (2*D)
        - 2D: dt < dx² / (4*D)
        - 3D: dt < dx² / (6*D)

    For implicit methods:
        Unconditionally stable, but accuracy degrades for large dt.
        Recommended: CFL < 10 for acceptable accuracy.

    **Physical Interpretation:**

    CFL = D*dt/dx² represents the ratio of:
        - Diffusion distance in time dt: sqrt(2*D*dt)
        - Mesh spacing: dx

    CFL << 1: Diffusion much slower than mesh resolution (well-resolved)
    CFL ~ 1:  Diffusion comparable to mesh resolution
    CFL >> 1: Diffusion faster than mesh resolution (poorly resolved)

    Examples
    --------
    Check stability for explicit heat diffusion (thermal diffusivity = 1e-6 m²/s):

    >>> is_stable, cfl, cfl_max = check_cfl_condition(
    ...     dt=0.1, dx=0.01, diffusivity=1e-6, method='explicit', ndim=2
    ... )
    >>> print(f"CFL = {cfl:.4f}, max = {cfl_max:.4f}, stable = {is_stable}")
    CFL = 0.0010, max = 0.0250, stable = True

    Check for implicit method (always stable but may have accuracy issues):

    >>> is_stable, cfl, cfl_max = check_cfl_condition(
    ...     dt=10, dx=0.01, diffusivity=1e-6, method='implicit'
    ... )
    >>> # Will warn if CFL > 10 (accuracy concern)

    See Also
    --------
    validate_mesh_quality : Check mesh quality metrics
    check_convergence : Monitor iterative solver convergence
    """
    # Compute CFL number
    cfl_number = diffusivity * dt / (dx**2)

    if method == 'explicit':
        # Explicit methods have stability limit
        cfl_max = 1.0 / (2.0 * ndim)
        is_stable = cfl_number <= cfl_max

        if not is_stable:
            warnings.warn(
                f"\n"
                f"╔════════════════════════════════════════════════════════════╗\n"
                f"║ CFL STABILITY CONDITION VIOLATED                           ║\n"
                f"╠════════════════════════════════════════════════════════════╣\n"
                f"║ Method:      {method:20s}                        ║\n"
                f"║ Dimensions:  {ndim}D                                           ║\n"
                f"║                                                            ║\n"
                f"║ CFL number:  {cfl_number:.6e} (current)                    ║\n"
                f"║ CFL limit:   {cfl_max:.6e} (maximum for stability)         ║\n"
                f"║                                                            ║\n"
                f"║ Time step:   dt = {dt:.6e} s                           ║\n"
                f"║ Mesh size:   dx = {dx:.6e} m                           ║\n"
                f"║ Diffusivity: D  = {diffusivity:.6e} m²/s                  ║\n"
                f"║                                                            ║\n"
                f"║ SOLUTION MAY BE UNSTABLE!                                  ║\n"
                f"║                                                            ║\n"
                f"║ Recommended maximum time step:                             ║\n"
                f"║   dt_max = {cfl_max * dx**2 / diffusivity:.6e} s           ║\n"
                f"╚════════════════════════════════════════════════════════════╝\n",
                UserWarning,
                stacklevel=2
            )

    elif method == 'implicit':
        # Implicit methods are unconditionally stable
        cfl_max = np.inf
        is_stable = True

        # But warn about accuracy if CFL is very large
        if cfl_number > 10:
            warnings.warn(
                f"\n"
                f"╔════════════════════════════════════════════════════════════╗\n"
                f"║ LARGE TIME STEP WARNING                                    ║\n"
                f"╠════════════════════════════════════════════════════════════╣\n"
                f"║ Method:      {method:20s} (unconditionally stable)   ║\n"
                f"║                                                            ║\n"
                f"║ CFL number:  {cfl_number:.6e}                              ║\n"
                f"║ Recommended: CFL < 10 for good accuracy                    ║\n"
                f"║                                                            ║\n"
                f"║ Time step:   dt = {dt:.6e} s                           ║\n"
                f"║ Mesh size:   dx = {dx:.6e} m                           ║\n"
                f"║                                                            ║\n"
                f"║ Solution is STABLE but may be INACCURATE.                  ║\n"
                f"║ Consider reducing time step for better accuracy.           ║\n"
                f"╚════════════════════════════════════════════════════════════╝\n",
                UserWarning,
                stacklevel=2
            )

    else:
        raise ValueError(f"Unknown method: '{method}'. Use 'explicit' or 'implicit'.")

    return is_stable, cfl_number, cfl_max


def validate_mesh_quality(GCOORD, EL2NOD):
    """
    Check triangular mesh quality metrics.

    Quality metrics help identify poorly-shaped elements that can degrade
    FEM accuracy and cause numerical issues.

    Parameters
    ----------
    GCOORD : ndarray, shape (nnod, 2)
        Node coordinates [m]
        GCOORD[i, 0] = x-coordinate of node i
        GCOORD[i, 1] = y-coordinate of node i
    EL2NOD : ndarray, shape (nel, 3)
        Element-to-node connectivity
        EL2NOD[iel, :] = [node1, node2, node3] for element iel

    Returns
    -------
    metrics : dict
        Dictionary containing mesh quality metrics:
            'min_angle' : float
                Minimum angle in mesh [degrees]
            'max_angle' : float
                Maximum angle in mesh [degrees]
            'min_quality' : float
                Minimum element quality [0-1]
            'mean_quality' : float
                Mean element quality [0-1]
            'aspect_ratios' : ndarray
                Aspect ratio for each element
            'poor_elements' : ndarray
                Indices of elements with quality < 0.3

    Notes
    -----
    **Element Quality Metric:**

    Quality Q is defined as:

        Q = 4 * sqrt(3) * A / (l1² + l2² + l3²)

    where:
        - A is the element area
        - l1, l2, l3 are the edge lengths

    Quality ranges from 0 to 1:
        - Q = 1.0: Equilateral triangle (ideal)
        - Q = 0.5: Right triangle with legs 1:1
        - Q < 0.3: Poor quality (avoid if possible)
        - Q → 0:   Degenerate (needle-like or flat)

    **Angle Constraints:**

    For good FEM performance:
        - Minimum angle > 20° (prefer > 30°)
        - Maximum angle < 120° (prefer < 90°)

    **Aspect Ratio:**

    Ratio of longest to shortest edge.
        - Ideal: ~1 (equilateral)
        - Acceptable: < 3
        - Poor: > 5

    Examples
    --------
    Check mesh quality:

    >>> import triangle as tr
    >>> # Create simple mesh
    >>> A = dict(vertices=[[0,0], [1,0], [1,1], [0,1]],
    ...          segments=[(0,1), (1,2), (2,3), (3,0)])
    >>> B = tr.triangulate(A, 'pq30a0.01')
    >>> GCOORD = B['vertices']
    >>> EL2NOD = B['triangles']
    >>> metrics = validate_mesh_quality(GCOORD, EL2NOD)
    >>> print(f"Min angle: {metrics['min_angle']:.1f}°")
    >>> print(f"Min quality: {metrics['min_quality']:.3f}")
    >>> if len(metrics['poor_elements']) > 0:
    ...     print(f"Warning: {len(metrics['poor_elements'])} poor elements")

    See Also
    --------
    mesh_utils.make_box : Create rectangular mesh
    mesh_utils.make_inclusion : Add circular inclusions
    """
    nel = EL2NOD.shape[0]

    # Initialize arrays
    qualities = np.zeros(nel)
    min_angles = np.zeros(nel)
    max_angles = np.zeros(nel)
    aspect_ratios = np.zeros(nel)

    for iel in range(nel):
        # Get element node coordinates
        nodes = EL2NOD[iel, :]
        coords = GCOORD[nodes, :]  # shape (3, 2)

        # Compute edge vectors
        v1 = coords[1, :] - coords[0, :]  # Edge 0->1
        v2 = coords[2, :] - coords[1, :]  # Edge 1->2
        v3 = coords[0, :] - coords[2, :]  # Edge 2->0

        # Edge lengths
        l1 = np.linalg.norm(v1)
        l2 = np.linalg.norm(v2)
        l3 = np.linalg.norm(v3)

        # Element area (using cross product)
        area = 0.5 * abs(np.cross(v1, -v3))

        # Quality metric
        perimeter_sq = l1**2 + l2**2 + l3**2
        qualities[iel] = 4 * np.sqrt(3) * area / perimeter_sq

        # Aspect ratio
        aspect_ratios[iel] = max(l1, l2, l3) / min(l1, l2, l3)

        # Angles (using law of cosines)
        # Angle at vertex 0
        angle0 = np.arccos(np.clip((l1**2 + l3**2 - l2**2) / (2*l1*l3), -1, 1))
        # Angle at vertex 1
        angle1 = np.arccos(np.clip((l1**2 + l2**2 - l3**2) / (2*l1*l2), -1, 1))
        # Angle at vertex 2
        angle2 = np.arccos(np.clip((l2**2 + l3**2 - l1**2) / (2*l2*l3), -1, 1))

        angles = np.array([angle0, angle1, angle2])
        min_angles[iel] = np.min(angles)
        max_angles[iel] = np.max(angles)

    # Convert angles to degrees
    min_angle_deg = np.rad2deg(np.min(min_angles))
    max_angle_deg = np.rad2deg(np.max(max_angles))

    # Find poor quality elements (Q < 0.3)
    poor_elements = np.where(qualities < 0.3)[0]

    # Package results
    metrics = {
        'min_angle': min_angle_deg,
        'max_angle': max_angle_deg,
        'min_quality': np.min(qualities),
        'mean_quality': np.mean(qualities),
        'aspect_ratios': aspect_ratios,
        'poor_elements': poor_elements,
        'qualities': qualities  # Full array for visualization
    }

    # Print warnings if mesh quality is poor
    if min_angle_deg < 20:
        warnings.warn(
            f"Mesh has small angles: min angle = {min_angle_deg:.1f}°. "
            f"Recommend min angle > 20° (preferably > 30°).",
            UserWarning
        )

    if max_angle_deg > 120:
        warnings.warn(
            f"Mesh has large angles: max angle = {max_angle_deg:.1f}°. "
            f"Recommend max angle < 120° (preferably < 90°).",
            UserWarning
        )

    if len(poor_elements) > 0:
        warnings.warn(
            f"Mesh has {len(poor_elements)} poor-quality elements (Q < 0.3). "
            f"This may affect solution accuracy.",
            UserWarning
        )

    return metrics


def check_convergence(residual, tolerance=1e-6, iteration=0, max_iter=100):
    """
    Check if iterative solver has converged.

    This function monitors convergence of nonlinear or iterative solvers,
    providing clear feedback on convergence status.

    Parameters
    ----------
    residual : float
        Current residual norm (measure of error)
    tolerance : float, optional
        Convergence tolerance. Default is 1e-6.
    iteration : int, optional
        Current iteration number. Default is 0.
    max_iter : int, optional
        Maximum allowed iterations. Default is 100.

    Returns
    -------
    converged : bool
        True if converged or diverged (stopping condition met)
    status : str
        Convergence status message:
            'converged' - Residual below tolerance
            'max_iter'  - Maximum iterations reached
            'diverged'  - Solution diverged (NaN or Inf)
            'iterating' - Still iterating

    Examples
    --------
    Monitor Newton iteration for nonlinear problem:

    >>> for it in range(max_iter):
    ...     # Solve one iteration
    ...     residual = compute_residual()
    ...
    ...     # Check convergence
    ...     converged, status = check_convergence(
    ...         residual, tolerance=1e-6, iteration=it, max_iter=20
    ...     )
    ...
    ...     print(f"Iter {it}: residual = {residual:.2e}, status = {status}")
    ...
    ...     if converged:
    ...         if status == 'converged':
    ...             print("✓ Solution converged!")
    ...         elif status == 'max_iter':
    ...             print("✗ Maximum iterations reached without convergence")
    ...         elif status == 'diverged':
    ...             print("✗ Solution diverged!")
    ...         break

    See Also
    --------
    check_cfl_condition : Check numerical stability
    """
    # Check for NaN or Inf (divergence)
    if np.isnan(residual) or np.isinf(residual):
        return True, 'diverged'

    # Check for convergence
    if residual < tolerance:
        return True, 'converged'

    # Check for maximum iterations
    if iteration >= max_iter - 1:
        return True, 'max_iter'

    # Still iterating
    return False, 'iterating'


def validate_parameters(params):
    """
    Validate model parameters before simulation.

    Checks that all required parameters are present and have physically
    reasonable values.

    Parameters
    ----------
    params : dict
        Dictionary of model parameters. Expected keys depend on problem type.
        Common keys:
            'dt' : float - Time step [s]
            'nt' : int - Number of time steps
            'diffusivity' or 'D' : float - Diffusion coefficient [m²/s]
            'conductivity' or 'k' : float - Thermal conductivity [W/(m·K)]

    Returns
    -------
    None
        Raises ValueError if parameters are invalid.
        Issues warnings for suspicious values.

    Raises
    ------
    ValueError
        If any parameter is invalid or out of physical range.

    Examples
    --------
    >>> params = {
    ...     'dt': 0.001,
    ...     'nt': 1000,
    ...     'D_A': 1.0,
    ...     'D_B': 0.5,
    ...     'f': 0.055,
    ...     'k': 0.062
    ... }
    >>> validate_parameters(params)  # Will pass or raise ValueError

    See Also
    --------
    check_cfl_condition : Check time step stability
    """
    # Check time stepping parameters
    if 'dt' in params:
        if params['dt'] <= 0:
            raise ValueError(f"Time step must be positive, got dt={params['dt']}")
        if params['dt'] > 1e6:
            warnings.warn(f"Very large time step: dt={params['dt']:.2e} s", UserWarning)

    if 'nt' in params:
        if params['nt'] <= 0:
            raise ValueError(f"Number of time steps must be positive, got nt={params['nt']}")
        if params['nt'] > 1e6:
            warnings.warn(f"Very large number of time steps: nt={params['nt']}", UserWarning)

    # Check diffusivity parameters
    for key in ['diffusivity', 'D', 'D_A', 'D_B', 'kappa']:
        if key in params:
            if params[key] < 0:
                raise ValueError(f"{key} must be non-negative, got {params[key]}")
            if params[key] == 0:
                warnings.warn(f"{key} = 0 (no diffusion)", UserWarning)

    # Check conductivity
    for key in ['conductivity', 'k', 'thermal_conductivity']:
        if key in params:
            if params[key] <= 0:
                raise ValueError(f"{key} must be positive, got {params[key]}")

    # Check reaction coefficients (if present)
    for key in ['f', 'feed_rate']:
        if key in params:
            if not (0 < params[key] < 1):
                warnings.warn(
                    f"{key}={params[key]} outside typical range [0, 1] for Gray-Scott model",
                    UserWarning
                )

    for key in ['k_coeff', 'kill_rate']:
        if key in params:
            if not (0 < params[key] < 1):
                warnings.warn(
                    f"{key}={params[key]} outside typical range [0, 1] for Gray-Scott model",
                    UserWarning
                )


if __name__ == "__main__":
    """
    Demonstration of validation utilities.
    """
    print("=" * 70)
    print("Validation Utilities - Demonstration")
    print("=" * 70)
    print()

    # Test 1: CFL condition check
    print("Test 1: CFL Stability Check")
    print("-" * 70)

    # Stable case
    print("Case 1a: Explicit method, stable")
    is_stable, cfl, cfl_max = check_cfl_condition(
        dt=0.01, dx=0.1, diffusivity=1e-6, method='explicit', ndim=2
    )
    print(f"  CFL = {cfl:.6f}, max = {cfl_max:.6f}, stable = {is_stable}")
    print()

    # Unstable case (will trigger warning)
    print("Case 1b: Explicit method, unstable (will show warning)")
    is_stable, cfl, cfl_max = check_cfl_condition(
        dt=100, dx=0.1, diffusivity=1e-6, method='explicit', ndim=2
    )
    print(f"  CFL = {cfl:.6f}, max = {cfl_max:.6f}, stable = {is_stable}")
    print()

    # Test 2: Parameter validation
    print("Test 2: Parameter Validation")
    print("-" * 70)

    params_good = {
        'dt': 0.001,
        'nt': 1000,
        'D_A': 1.0,
        'D_B': 0.5
    }
    print("Validating good parameters...")
    validate_parameters(params_good)
    print("  ✓ All parameters valid")
    print()

    # Test 3: Convergence monitoring
    print("Test 3: Convergence Monitoring")
    print("-" * 70)
    residuals = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
    for it, res in enumerate(residuals):
        converged, status = check_convergence(res, tolerance=1e-6, iteration=it)
        print(f"  Iter {it}: residual = {res:.2e}, status = {status:10s}")
        if converged:
            break
    print()

    print("=" * 70)
    print("All demonstrations complete!")
    print("=" * 70)
