"""
Shape Functions for Triangular Finite Elements

This module provides shape functions and their derivatives for triangular
elements used in the Finite Element Method (FEM). It supports both linear
(3-node) and quadratic (6-node) triangular elements.

The shape functions are defined in the reference element with vertices at:
    Node 1: (xi=0, eta=0)
    Node 2: (xi=1, eta=0)
    Node 3: (xi=0, eta=1)

For 6-node elements, additional mid-side nodes are located at:
    Node 4: Midpoint of edge 2-3
    Node 5: Midpoint of edge 1-3
    Node 6: Midpoint of edge 1-2

Mathematical Background
-----------------------
The shape functions satisfy:
    1. Partition of unity: sum(N_i) = 1 for all (xi, eta)
    2. Kronecker delta property: N_i(node_j) = delta_ij
    3. Completeness: Can represent polynomials up to a certain degree

For linear elements (nnodel=3):
    N_i are degree-1 polynomials (can represent linear functions)

For quadratic elements (nnodel=6):
    N_i are degree-2 polynomials (can represent quadratic functions)

References
----------
.. [1] Zienkiewicz, O. C., & Taylor, R. L. (2000).
       The Finite Element Method, Volume 1: The Basis. 5th ed.
.. [2] Hughes, T. J. R. (2000).
       The Finite Element Method: Linear Static and Dynamic Finite Element Analysis.

Author: HTF Lecture Series
Date: 2025
"""

import numpy as np


def shapes_tri(xi, eta, nnodel=3):
    """
    Compute shape functions and derivatives for triangular elements.

    This function evaluates the shape functions and their derivatives
    with respect to local coordinates (xi, eta) for triangular elements.

    Parameters
    ----------
    xi : float
        First local coordinate (0 <= xi <= 1)
    eta : float
        Second local coordinate (0 <= eta <= 1)
        Constraint: xi + eta <= 1 (inside reference triangle)
    nnodel : int, optional
        Number of nodes per element. Default is 3.
        Supported values:
            - 3: Linear triangular element (P1)
            - 6: Quadratic triangular element (P2)

    Returns
    -------
    N : ndarray, shape (nnodel,)
        Shape function values at the specified point.
        N[i] is the value of shape function i at (xi, eta)
    dNds : ndarray, shape (2, nnodel)
        Shape function derivatives with respect to local coordinates.
        dNds[0, i] = dN_i/dxi   (derivative w.r.t. first local coordinate)
        dNds[1, i] = dN_i/deta  (derivative w.r.t. second local coordinate)

    Notes
    -----
    The local coordinates (xi, eta) use area coordinates (also called
    barycentric coordinates):
        eta1 = 1 - xi - eta  (associated with node 1)
        eta2 = xi            (associated with node 2)
        eta3 = eta           (associated with node 3)

    These coordinates represent the normalized distance from each edge:
        eta1 + eta2 + eta3 = 1 (partition of unity)

    Linear Element Shape Functions (nnodel=3):
        N1 = 1 - xi - eta
        N2 = xi
        N3 = eta

    Quadratic Element Shape Functions (nnodel=6):
        N1 = (1-xi-eta)*(2*(1-xi-eta)-1)     Corner node 1
        N2 = xi*(2*xi-1)                      Corner node 2
        N3 = eta*(2*eta-1)                    Corner node 3
        N4 = 4*xi*eta                         Mid-side node 2-3
        N5 = 4*(1-xi-eta)*eta                 Mid-side node 1-3
        N6 = 4*(1-xi-eta)*xi                  Mid-side node 1-2

    Raises
    ------
    ValueError
        If nnodel is not 3 or 6

    Examples
    --------
    Evaluate shape functions at the centroid of the triangle:

    >>> N, dNds = shapes_tri(1/3, 1/3, nnodel=3)
    >>> print(N)
    [0.33333333 0.33333333 0.33333333]
    >>> print(dNds)
    [[-1.  1.  0.]
     [-1.  0.  1.]]

    For a 6-node element at the same point:

    >>> N, dNds = shapes_tri(1/3, 1/3, nnodel=6)
    >>> print(np.sum(N))  # Should be 1 (partition of unity)
    1.0

    Evaluate at node 1 (corner):

    >>> N, dNds = shapes_tri(0, 0, nnodel=3)
    >>> print(N)  # N1=1, others=0 (Kronecker delta property)
    [1. 0. 0.]

    See Also
    --------
    fem_utils.compute_jacobian : Compute Jacobian matrix using dNds
    fem_integration.get_triangle_quadrature : Get integration points for element assembly
    """

    # Define area/barycentric coordinates
    eta2 = xi
    eta3 = eta
    eta1 = 1 - eta2 - eta3

    if nnodel == 3:
        # ========================================
        # Linear Triangular Element (P1)
        # ========================================
        # 3 corner nodes

        # Shape functions (degree 1 polynomials)
        N1 = eta1
        N2 = eta2
        N3 = eta3

        N = np.array([N1, N2, N3])

        # Shape function derivatives (constant for linear elements)
        dNds = np.zeros((2, 3))

        # Derivatives of N1 = 1 - xi - eta
        dNds[0, 0] = -1  # dN1/dxi
        dNds[1, 0] = -1  # dN1/deta

        # Derivatives of N2 = xi
        dNds[0, 1] = 1   # dN2/dxi
        dNds[1, 1] = 0   # dN2/deta

        # Derivatives of N3 = eta
        dNds[0, 2] = 0   # dN3/dxi
        dNds[1, 2] = 1   # dN3/deta

    elif nnodel == 6:
        # ========================================
        # Quadratic Triangular Element (P2)
        # ========================================
        # 3 corner nodes + 3 mid-side nodes

        # Shape functions (degree 2 polynomials)
        # Corner nodes (include self-interaction term)
        N1 = eta1 * (2*eta1 - 1)
        N2 = eta2 * (2*eta2 - 1)
        N3 = eta3 * (2*eta3 - 1)

        # Mid-side nodes (bilinear interaction terms)
        N4 = 4 * eta2 * eta3  # Mid-side node between nodes 2-3
        N5 = 4 * eta1 * eta3  # Mid-side node between nodes 1-3
        N6 = 4 * eta1 * eta2  # Mid-side node between nodes 1-2

        N = np.array([N1, N2, N3, N4, N5, N6])

        # Shape function derivatives
        dNds = np.zeros((2, 6))

        # Derivatives of corner node shape functions
        # N1 = eta1*(2*eta1-1), where eta1 = 1-xi-eta
        # dN1/dxi = d(eta1*(2*eta1-1))/dxi = (2*eta1-1)*(-1) + eta1*2*(-1) = -4*eta1+1
        dNds[0, 0] = 1 - 4*eta1  # dN1/dxi
        dNds[1, 0] = 1 - 4*eta1  # dN1/deta

        # N2 = eta2*(2*eta2-1), where eta2 = xi
        # dN2/dxi = 4*xi - 1
        dNds[0, 1] = -1 + 4*eta2  # dN2/dxi = 4*xi - 1
        dNds[1, 1] = 0            # dN2/deta = 0

        # N3 = eta3*(2*eta3-1), where eta3 = eta
        # dN3/deta = 4*eta - 1
        dNds[0, 2] = 0            # dN3/dxi = 0
        dNds[1, 2] = -1 + 4*eta3  # dN3/deta = 4*eta - 1

        # Derivatives of mid-side node shape functions
        # N4 = 4*xi*eta
        dNds[0, 3] = 4*eta3   # dN4/dxi = 4*eta
        dNds[1, 3] = 4*eta2   # dN4/deta = 4*xi

        # N5 = 4*(1-xi-eta)*eta = 4*eta1*eta3
        dNds[0, 4] = -4*eta3  # dN5/dxi = -4*eta
        dNds[1, 4] = 4*eta1 - 4*eta3  # dN5/deta = 4*(1-xi-eta) - 4*eta = 4 - 4*xi - 8*eta

        # N6 = 4*(1-xi-eta)*xi = 4*eta1*eta2
        dNds[0, 5] = 4*eta1 - 4*eta2  # dN6/dxi = 4*(1-xi-eta) - 4*xi = 4 - 8*xi - 4*eta
        dNds[1, 5] = -4*eta2  # dN6/deta = -4*xi

    else:
        raise ValueError(f"Unsupported number of nodes: nnodel={nnodel}. "
                         f"Supported values are 3 (linear) and 6 (quadratic).")

    return N, dNds


def verify_shape_functions(nnodel=3, tol=1e-10):
    """
    Verify mathematical properties of shape functions.

    This utility function checks that the shape functions satisfy
    fundamental properties required for finite element analysis.

    Parameters
    ----------
    nnodel : int, optional
        Number of nodes per element (3 or 6). Default is 3.
    tol : float, optional
        Tolerance for numerical comparisons. Default is 1e-10.

    Returns
    -------
    bool
        True if all tests pass, False otherwise.

    Tests Performed
    ---------------
    1. Partition of Unity: sum(N_i) = 1 at random points
    2. Kronecker Delta: N_i(node_j) = 1 if i==j, else 0
    3. Derivatives: Check analytical derivative values

    Examples
    --------
    >>> verify_shape_functions(nnodel=3)
    Testing shape functions for nnodel=3...
    ✓ Partition of unity satisfied
    ✓ Kronecker delta property satisfied
    ✓ All tests passed!
    True
    """
    print(f"Testing shape functions for nnodel={nnodel}...")

    # Test 1: Partition of unity at random points
    np.random.seed(42)
    for _ in range(10):
        xi = np.random.uniform(0, 1)
        eta = np.random.uniform(0, 1-xi)  # Ensure inside triangle

        N, dNds = shapes_tri(xi, eta, nnodel)

        if abs(np.sum(N) - 1.0) > tol:
            print(f"✗ Partition of unity failed at (xi={xi:.3f}, eta={eta:.3f}): sum(N)={np.sum(N):.10f}")
            return False

    print("✓ Partition of unity satisfied")

    # Test 2: Kronecker delta property at nodes
    if nnodel == 3:
        node_coords = [(0, 0), (1, 0), (0, 1)]
    elif nnodel == 6:
        node_coords = [(0, 0), (1, 0), (0, 1), (0.5, 0.5), (0, 0.5), (0.5, 0)]
    else:
        print(f"✗ Unsupported nnodel={nnodel}")
        return False

    for i, (xi, eta) in enumerate(node_coords):
        N, dNds = shapes_tri(xi, eta, nnodel)

        # N[i] should be 1, all others should be 0
        expected = np.zeros(nnodel)
        expected[i] = 1.0

        if np.linalg.norm(N - expected) > tol:
            print(f"✗ Kronecker delta failed at node {i}: N={N}")
            return False

    print("✓ Kronecker delta property satisfied")

    # Test 3: Derivatives for linear elements should be constant
    if nnodel == 3:
        N1, dNds1 = shapes_tri(0.2, 0.3, nnodel)
        N2, dNds2 = shapes_tri(0.5, 0.1, nnodel)

        if np.linalg.norm(dNds1 - dNds2) > tol:
            print("✗ Linear element derivatives are not constant")
            return False

        print("✓ Linear element derivatives are constant")

    print("✓ All tests passed!")
    return True


if __name__ == "__main__":
    # Run verification tests when module is executed directly
    print("=" * 60)
    print("FEM Shape Functions - Verification Tests")
    print("=" * 60)
    print()

    verify_shape_functions(nnodel=3)
    print()
    verify_shape_functions(nnodel=6)
    print()

    # Example usage
    print("=" * 60)
    print("Example Usage")
    print("=" * 60)
    print()
    print("Evaluate at triangle centroid (xi=1/3, eta=1/3):")
    N, dNds = shapes_tri(1/3, 1/3, nnodel=3)
    print(f"Shape functions N = {N}")
    print(f"Derivatives dNds =\n{dNds}")
