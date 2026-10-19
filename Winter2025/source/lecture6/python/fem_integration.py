"""
Gauss Quadrature Rules for Triangular Elements

This module provides integration points and weights for numerical integration
over triangular elements using Gauss quadrature rules. These rules are essential
for finite element method (FEM) calculations where integrals must be evaluated
numerically.

Quadrature rules are specified by their order (number of integration points),
which determines the polynomial degree that can be integrated exactly.

Reference
---------
Based on lecture7/python/ip_triangle.py

Author: HTF Lecture Series
Date: 2025
"""

import numpy as np


def ip_triangle(nip):
    """
    Get Gauss quadrature points and weights for triangular elements.

    This function returns integration points and weights for numerical
    integration over the reference triangle with vertices at (0,0), (1,0),
    and (0,1).

    Parameters
    ----------
    nip : int
        Number of integration points.
        Supported values: 1, 3, 4, 6, 7, 12, 25, 77 (Labatto)

    Returns
    -------
    ipx : ndarray, shape (nip, 2)
        Integration point coordinates in reference element.
        ipx[i, 0] = xi coordinate of point i
        ipx[i, 1] = eta coordinate of point i
        Each point satisfies: 0 <= xi, eta and xi + eta <= 1
    ipw : ndarray, shape (nip,)
        Integration weights.
        Sum of weights equals 0.5 (area of reference triangle)

    Notes
    -----
    **Integration Formula:**

    For a function f(xi, eta) over the reference triangle:

        ∫∫_T f(xi, eta) dA ≈ sum_{i=1}^{nip} w_i * f(xi_i, eta_i) * detJ

    where w_i are weights, (xi_i, eta_i) are integration points,
    and detJ is the Jacobian determinant.

    **Polynomial Exactness:**

    Each quadrature rule integrates polynomials up to a certain degree exactly:

    - nip=1:  degree 1 (linear functions)
    - nip=3:  degree 2 (quadratic functions) - most common
    - nip=4:  degree 3 (cubic functions)
    - nip=6:  degree 4
    - nip=7:  degree 5
    - nip=12: degree 6
    - nip=25: degree 10
    - nip=77: Labatto rule (includes vertices and edge midpoints)

    **Usage Guidelines:**

    - Linear elements (P1): Use nip=3 (exact for degree 2)
    - Quadratic elements (P2): Use nip=6 or 7 (exact for degree 4-5)
    - Higher-order elements: Use nip=12 or 25

    **Reference Triangle:**

    The integration points are defined on the reference triangle:
        Vertex 1: (0, 0)
        Vertex 2: (1, 0)
        Vertex 3: (0, 1)

    which has area = 0.5, hence sum(weights) = 0.5.

    Raises
    ------
    ValueError
        If nip is not one of the supported values.

    Examples
    --------
    Get 3-point quadrature rule (most common):

    >>> ipx, ipw = ip_triangle(nip=3)
    >>> print(f"Integration points:\n{ipx}")
    [[0.16666667 0.16666667]
     [0.66666667 0.16666667]
     [0.16666667 0.66666667]]
    >>> print(f"Weights: {ipw}")
    [0.16666667 0.16666667 0.16666667]
    >>> print(f"Sum of weights: {np.sum(ipw):.6f}")
    Sum of weights: 0.500000

    Verify integration of constant function (should equal area = 0.5):

    >>> ipx, ipw = ip_triangle(3)
    >>> integral = np.sum(ipw)  # f(x,y) = 1
    >>> print(f"Integral of f=1: {integral:.6f}")
    Integral of f=1: 0.500000

    Verify integration of linear function f(xi, eta) = xi:

    >>> ipx, ipw = ip_triangle(3)
    >>> integral = np.sum(ipw * ipx[:, 0])
    >>> exact = 1/6  # int_0^1 int_0^(1-xi) xi deta dxi
    >>> print(f"Numerical: {integral:.6f}, Exact: {exact:.6f}")
    Numerical: 0.166667, Exact: 0.166667

    Use in FEM element assembly:

    >>> from fem_shapes import shapes_tri
    >>> ipx, ipw = ip_triangle(nip=3)
    >>> for ip in range(len(ipw)):
    ...     xi, eta = ipx[ip, 0], ipx[ip, 1]
    ...     N, dNds = shapes_tri(xi, eta, nnodel=3)
    ...     # ... perform integration with weight ipw[ip]

    See Also
    --------
    fem_shapes.shapes_tri : Evaluate shape functions at integration points
    fem_utils.compute_jacobian : Compute Jacobian for integration
    """

    if nip == 1:
        # =====================================================
        # 1-point rule (degree 1 precision)
        # =====================================================
        # Single point at triangle centroid
        ipx = np.array([[1/3, 1/3]])
        ipw = np.array([0.5])

    elif nip == 3:
        # =====================================================
        # 3-point rule (degree 2 precision) - MOST COMMON
        # =====================================================
        # Points at 1/6 and 2/3 along triangle edges
        # This is the standard rule for linear elements
        ipx = np.array([
            [1/6, 1/6],   # Near vertex 1
            [2/3, 1/6],   # Near vertex 2
            [1/6, 2/3]    # Near vertex 3
        ])
        ipw = np.array([1/6, 1/6, 1/6])

    elif nip == 4:
        # =====================================================
        # 4-point rule (degree 3 precision)
        # =====================================================
        # Includes centroid (negative weight) + 3 points
        ipx = np.array([
            [1.0/3.0, 1.0/3.0],  # Centroid
            [0.6, 0.2],
            [0.2, 0.6],
            [0.2, 0.2]
        ])
        ipw = np.array([-27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0]) * 0.5

    elif nip == 6:
        # =====================================================
        # 6-point rule (degree 4 precision)
        # =====================================================
        # Two groups of 3 points each (3-fold symmetry)
        g1 = (8 - np.sqrt(10) + np.sqrt(38 - 44 * np.sqrt(2/5))) / 18
        g2 = (8 - np.sqrt(10) - np.sqrt(38 - 44 * np.sqrt(2/5))) / 18

        ipx = np.array([
            [1 - 2 * g1, g1],
            [g1, 1 - 2 * g1],
            [g1, g1],
            [1 - 2 * g2, g2],
            [g2, 1 - 2 * g2],
            [g2, g2]
        ])

        w1 = (620 + np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
        w2 = (620 - np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
        ipw = np.array([w1, w1, w1, w2, w2, w2]) * 0.5

    elif nip == 7:
        # =====================================================
        # 7-point rule (degree 5 precision) - RECOMMENDED FOR P2
        # =====================================================
        # Centroid + two groups of 3 points
        g1 = (6 - np.sqrt(15)) / 21
        g2 = (6 + np.sqrt(15)) / 21

        ipx = np.array([
            [1/3, 1/3],           # Centroid
            [1 - 2 * g1, g1],     # Group 1
            [g1, 1 - 2 * g1],
            [g1, g1],
            [1 - 2 * g2, g2],     # Group 2
            [g2, 1 - 2 * g2],
            [g2, g2]
        ])

        w1 = (155 - np.sqrt(15)) / 1200
        w2 = (155 + np.sqrt(15)) / 1200
        ipw = np.array([0.225, w1, w1, w1, w2, w2, w2]) * 0.5

    elif nip == 77:
        # =====================================================
        # 7-point Labatto rule (includes vertices)
        # =====================================================
        # Special rule with integration points at vertices,
        # edge midpoints, and centroid
        ipx = np.array([
            [0, 0],       # Vertex 1
            [1, 0],       # Vertex 2
            [0, 1],       # Vertex 3
            [1/2, 1/2],   # Edge midpoint 2-3
            [1/2, 0],     # Edge midpoint 1-2
            [0, 1/2],     # Edge midpoint 1-3
            [1/3, 1/3]    # Centroid
        ])
        ipw = np.array([1/40, 1/40, 1/40, 1/15, 1/15, 1/15, 9/40])

    elif nip == 12:
        # =====================================================
        # 12-point rule (degree 6 precision)
        # =====================================================
        # Three groups of 3 or 4 points with 3-fold symmetry
        ipx = np.array([
            [0.873821971016996, 0.063089014491502],
            [0.063089014491502, 0.873821971016996],
            [0.063089014491502, 0.063089014491502],
            [0.501426509658179, 0.249286745170910],
            [0.249286745170910, 0.501426509658179],
            [0.249286745170910, 0.249286745170910],
            [0.636502499121399, 0.310352451033785],
            [0.310352451033785, 0.636502499121399],
            [0.636502499121399, 0.053145049844816],
            [0.053145049844816, 0.636502499121399],
            [0.310352451033785, 0.053145049844816],
            [0.053145049844816, 0.310352451033785]
        ])

        ipw = np.array([
            0.050844906370207, 0.050844906370207, 0.050844906370207,
            0.116786275726379, 0.116786275726379, 0.116786275726379,
            0.082851075618374, 0.082851075618374, 0.082851075618374,
            0.082851075618374, 0.082851075618374, 0.082851075618374
        ]) * 0.5

    elif nip == 25:
        # =====================================================
        # 25-point rule (degree 10 precision) - HIGH ACCURACY
        # =====================================================
        # Multiple groups of points for very high order integration
        ipx = np.array([
            [0.333333333333333, 0.333333333333333],
            [0.028844733232685, 0.485577633383657],
            [0.485577633383657, 0.028844733232685],
            [0.485577633383657, 0.485577633383657],
            [0.781036849029926, 0.109481575485037],
            [0.109481575485037, 0.781036849029926],
            [0.109481575485037, 0.109481575485037],
            [0.141707219414880, 0.307939838764121],
            [0.307939838764121, 0.141707219414880],
            [0.141707219414880, 0.550352941820999],
            [0.550352941820999, 0.141707219414880],
            [0.307939838764121, 0.550352941820999],
            [0.550352941820999, 0.307939838764121],
            [0.025003534762686, 0.246672560639903],
            [0.246672560639903, 0.025003534762686],
            [0.025003534762686, 0.728323904597411],
            [0.728323904597411, 0.025003534762686],
            [0.246672560639903, 0.728323904597411],
            [0.728323904597411, 0.246672560639903],
            [0.009540815400299, 0.066803251012200],
            [0.066803251012200, 0.009540815400299],
            [0.009540815400299, 0.923655933587500],
            [0.923655933587500, 0.009540815400299],
            [0.066803251012200, 0.923655933587500],
            [0.923655933587500, 0.066803251012200]
        ])

        ipw = np.array([
            0.090817990382754,
            0.036725957756467, 0.036725957756467, 0.036725957756467,
            0.045321059435528, 0.045321059435528, 0.045321059435528,
            0.072757916845420, 0.072757916845420, 0.072757916845420,
            0.072757916845420, 0.072757916845420, 0.072757916845420,
            0.028327242531057, 0.028327242531057, 0.028327242531057,
            0.028327242531057, 0.028327242531057, 0.028327242531057,
            0.009421666963733, 0.009421666963733, 0.009421666963733,
            0.009421666963733, 0.009421666963733, 0.009421666963733
        ]) * 0.5

    else:
        raise ValueError(
            f"Unsupported number of integration points: nip={nip}.\n"
            f"Supported values: 1, 3, 4, 6, 7, 12, 25, 77 (Labatto)"
        )

    return ipx, ipw


def int_points_triangle(nip):
    """
    Legacy wrapper for ip_triangle() for backward compatibility.

    This function provides the same functionality as ip_triangle() but
    returns values in a different format (transposed) to match the old
    int_points_triangle.py interface.

    Parameters
    ----------
    nip : int
        Number of integration points (3 or 6 supported)

    Returns
    -------
    gauss : ndarray, shape (nip, 2)
        Integration point coordinates (same as ip_triangle)
    weights : ndarray, shape (nip,)
        Integration weights (same as ip_triangle)

    Notes
    -----
    This is a legacy function maintained for backward compatibility with
    existing code. New code should use ip_triangle() directly.

    For nip=3, returns the same 3-point Gauss rule.
    For nip=6, returns the 6-point rule.

    Examples
    --------
    >>> gauss, weights = int_points_triangle(nip=3)
    >>> # Same as: ipx, ipw = ip_triangle(3)
    """
    if nip not in [3, 6]:
        print(f"WARNING: int_points_triangle only supports nip=3 or 6.")
        print(f"         For nip={nip}, use ip_triangle() instead.")
        raise ValueError(f"nip={nip} not supported in legacy function")

    return ip_triangle(nip)


if __name__ == "__main__":
    """
    Demonstration and testing of integration rules.
    """
    print("=" * 70)
    print("Triangle Integration Rules - Demonstration")
    print("=" * 70)
    print()

    # Test all supported quadrature rules
    nip_values = [1, 3, 4, 6, 7, 12, 25, 77]

    print("Testing all supported quadrature rules:")
    print("-" * 70)

    for nip in nip_values:
        ipx, ipw = ip_triangle(nip)

        # Check sum of weights (should be 0.5)
        weight_sum = np.sum(ipw)

        # Check that all points are inside reference triangle
        inside = np.all(ipx >= 0) and np.all(ipx[:, 0] + ipx[:, 1] <= 1)

        print(f"nip={nip:2d}: {len(ipx):2d} points, "
              f"sum(w)={weight_sum:.10f}, "
              f"inside={'✓' if inside else '✗'}")

    print()

    # Test integration accuracy
    print("Testing integration accuracy:")
    print("-" * 70)

    def test_polynomial(degree, nip):
        """Test integration of polynomial x^p * y^q where p+q=degree"""
        ipx, ipw = ip_triangle(nip)

        # Test f(x,y) = x^degree
        if degree == 0:
            # Constant function: f = 1
            numerical = np.sum(ipw)
            exact = 0.5  # Area of reference triangle
        elif degree == 1:
            # Linear: f = x
            numerical = np.sum(ipw * ipx[:, 0])
            exact = 1/6  # int_0^1 int_0^(1-x) x dy dx
        elif degree == 2:
            # Quadratic: f = x^2
            numerical = np.sum(ipw * ipx[:, 0]**2)
            exact = 1/12  # int_0^1 int_0^(1-x) x^2 dy dx
        else:
            # General: f = x^degree
            numerical = np.sum(ipw * ipx[:, 0]**degree)
            exact = 1 / ((degree + 1) * (degree + 2))

        error = abs(numerical - exact)
        return error

    test_cases = [
        (3, [0, 1, 2]),      # 3-point: exact for degree 0-2
        (6, [0, 1, 2, 3, 4]), # 6-point: exact for degree 0-4
        (7, [0, 1, 2, 3, 4, 5]), # 7-point: exact for degree 0-5
    ]

    for nip, degrees in test_cases:
        errors = [test_polynomial(d, nip) for d in degrees]
        max_error = max(errors)
        status = "✓ EXACT" if max_error < 1e-12 else f"✗ error={max_error:.2e}"
        print(f"nip={nip:2d}: degree {degrees[0]}-{degrees[-1]} {status}")

    print()
    print("=" * 70)
    print("All demonstrations complete!")
    print("=" * 70)
