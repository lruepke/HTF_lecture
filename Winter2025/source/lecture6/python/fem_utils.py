"""
FEM Utility Functions for Finite Element Assembly

This module provides utility functions for finite element method (FEM)
computations, focusing on assembly operations, Jacobian transformations,
and other common FEM tasks.

The utilities are designed to be educational - prioritizing clarity
and explicitness over conciseness.

Author: HTF Lecture Series
Date: 2025
"""

import numpy as np
from scipy.sparse import csr_matrix


def compute_jacobian(dNds, ECOORD):
    """
    Compute Jacobian matrix for isoparametric mapping.

    The Jacobian matrix relates derivatives in the reference element
    coordinates (xi, eta) to derivatives in physical element coordinates (x, y).

    Parameters
    ----------
    dNds : ndarray, shape (2, nnodel)
        Shape function derivatives with respect to local coordinates.
        dNds[0, i] = dN_i/dxi   (derivative w.r.t. first local coordinate)
        dNds[1, i] = dN_i/deta  (derivative w.r.t. second local coordinate)
    ECOORD : ndarray, shape (nnodel, 2)
        Physical coordinates of element nodes.
        ECOORD[i, 0] = x-coordinate of node i
        ECOORD[i, 1] = y-coordinate of node i

    Returns
    -------
    Jac : ndarray, shape (2, 2)
        Jacobian matrix.
        Jac[0, 0] = dx/dxi,  Jac[0, 1] = dy/dxi
        Jac[1, 0] = dx/deta, Jac[1, 1] = dy/deta
    invJ : ndarray, shape (2, 2)
        Inverse Jacobian matrix.
        Used to transform derivatives: dN/dx = invJ * dN/dxi
    detJ : float
        Jacobian determinant.
        For 2D elements: detJ = 2 * element_area
        Used as weighting factor in numerical integration.

    Notes
    -----
    The Jacobian matrix is computed as:

        Jac = dNds @ ECOORD

    where @ denotes matrix multiplication. Explicitly:

        [dx/dxi   dy/dxi ]     [dN1/dxi   dN2/dxi   dN3/dxi ]   [x1  y1]
        [dx/deta  dy/deta]  =  [dN1/deta  dN2/deta  dN3/deta] * [x2  y2]
                                                                [x3  y3]

    The determinant detJ represents twice the element area for triangular
    elements. For a valid (non-degenerate) element, detJ > 0.

    Warnings
    --------
    If detJ ≤ 0, the element is inverted or degenerate. This typically
    indicates a mesh quality issue.

    Examples
    --------
    For an equilateral triangle with unit sides:

    >>> from fem_shapes import shapes_tri
    >>> ECOORD = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3)/2]])
    >>> N, dNds = shapes_tri(1/3, 1/3, nnodel=3)
    >>> Jac, invJ, detJ = compute_jacobian(dNds, ECOORD)
    >>> print(f"Element area = {detJ/2:.4f}")
    Element area = 0.4330

    See Also
    --------
    fem_shapes.shapes_tri : Compute shape functions and derivatives
    """
    # Compute Jacobian matrix: J = dNds * ECOORD
    # Shape: [2, nnodel] * [nnodel, 2] = [2, 2]
    Jac = np.matmul(dNds, ECOORD)

    # Compute inverse Jacobian
    # Used to transform derivatives from local to global coordinates
    invJ = np.linalg.inv(Jac)

    # Compute determinant
    # For 2D triangular elements: detJ = 2 * area
    detJ = np.linalg.det(Jac)

    # Check for degenerate elements
    if detJ <= 0:
        import warnings
        warnings.warn(
            f"Negative or zero Jacobian determinant: detJ = {detJ:.6e}. "
            "Element may be inverted or degenerate. Check mesh quality.",
            RuntimeWarning
        )

    return Jac, invJ, detJ


def create_assembly_indices_2dof(EL2NOD, nnodel):
    """
    Create assembly index arrays for coupled 2-DOF-per-node systems.

    This function generates the row and column indices (I, J) for assembling
    element matrices into the global sparse matrix for systems where each
    node has 2 degrees of freedom (DOFs). Common in coupled problems like
    reaction-diffusion systems with species A and B.

    This replaces the complex nested operations found in scripts like
    double_diff.py:136-138 with a cleaner, more readable implementation.

    Parameters
    ----------
    EL2NOD : ndarray, shape (nel, nnodel)
        Element-to-node connectivity matrix.
        EL2NOD[iel, i] = global node number for local node i of element iel
    nnodel : int
        Number of nodes per element (typically 3 for triangular elements)

    Returns
    -------
    I : ndarray, shape (nel, 2*nnodel*nnodel)
        Row indices for sparse matrix assembly.
        For each element, contains indices for both DOF blocks (A and B)
    J : ndarray, shape (nel, 2*nnodel*nnodel)
        Column indices for sparse matrix assembly.
        Paired with I to define matrix entry locations

    Notes
    -----
    For a coupled system with 2 DOFs per node (e.g., species A and B):
        - Node i has DOFs: [A_i, B_i] at global indices [2*i, 2*i+1]
        - The global matrix has structure:
            [A-A block | A-B block]
            [B-A block | B-B block]

    This function creates indices for all four blocks simultaneously.

    The assembly indices allow converting element matrices into sparse COO format:
        A_global = csr_matrix((K.flat, (I.flat, J.flat)), shape=(sdof, sdof))

    Examples
    --------
    For a single triangular element connecting nodes [5, 7, 12]:

    >>> EL2NOD = np.array([[5, 7, 12]])
    >>> I, J = create_assembly_indices_2dof(EL2NOD, nnodel=3)
    >>> print("A-DOF indices:", I[0, :9])  # First block (A-A)
    A-DOF indices: [10 10 10 14 14 14 24 24 24]
    >>> print("B-DOF indices:", I[0, 9:])  # Second block (B-B)
    B-DOF indices: [11 11 11 15 15 15 25 25 25]

    Performance Comparison
    ----------------------
    OLD (complex nested operations):
        I = np.concatenate((
            np.outer(2*EL2NOD[iel,:], np.ones(nnodel, dtype=int)).reshape(nnodel*nnodel),
            np.outer(2*EL2NOD[iel,:]+1, np.ones(nnodel, dtype=int)).reshape(nnodel*nnodel)
        ))

    NEW (this function):
        I, J = create_assembly_indices_2dof(EL2NOD, nnodel)

    See Also
    --------
    scipy.sparse.csr_matrix : Sparse matrix construction using COO format
    assemble_sparse_matrix_2dof : Assemble global matrix from element contributions
    """
    nel = EL2NOD.shape[0]

    # Pre-allocate index arrays
    # Total size: 2 blocks (A and B) * nnodel * nnodel entries per block
    I = np.zeros((nel, 2*nnodel*nnodel), dtype=int)
    J = np.zeros((nel, 2*nnodel*nnodel), dtype=int)

    # Loop over elements
    for iel in range(nel):
        # Get global node numbers for this element
        nodes = EL2NOD[iel, :]  # shape: (nnodel,)

        # Compute global DOF indices for both species
        dof_A = 2 * nodes      # A DOFs: [2*n1, 2*n2, 2*n3]
        dof_B = 2 * nodes + 1  # B DOFs: [2*n1+1, 2*n2+1, 2*n3+1]

        # ==========================================
        # Block 1: A-A coupling (first nnodel^2 entries)
        # ==========================================
        # Row indices: each A-DOF repeated nnodel times
        # For nnodel=3: [A1, A1, A1, A2, A2, A2, A3, A3, A3]
        I_AA = np.repeat(dof_A, nnodel)

        # Column indices: all A-DOFs repeated nnodel times
        # For nnodel=3: [A1, A2, A3, A1, A2, A3, A1, A2, A3]
        J_AA = np.tile(dof_A, nnodel)

        # ==========================================
        # Block 2: B-B coupling (second nnodel^2 entries)
        # ==========================================
        I_BB = np.repeat(dof_B, nnodel)
        J_BB = np.tile(dof_B, nnodel)

        # Concatenate both blocks
        # Total size: 2 * nnodel * nnodel
        I[iel, :] = np.concatenate([I_AA, I_BB])
        J[iel, :] = np.concatenate([J_AA, J_BB])

    return I, J


def assemble_sparse_matrix_2dof(I, J, K, sdof):
    """
    Assemble global sparse matrix from element contributions (2 DOF per node).

    Converts element-wise contributions into a global sparse matrix using
    coordinate (COO) format, then converts to compressed sparse row (CSR)
    format for efficient linear algebra operations.

    Parameters
    ----------
    I : ndarray, shape (nel, 2*nnodel*nnodel)
        Row indices from create_assembly_indices_2dof
    J : ndarray, shape (nel, 2*nnodel*nnodel)
        Column indices from create_assembly_indices_2dof
    K : ndarray, shape (nel, 2*nnodel*nnodel)
        Element stiffness/mass matrix values.
        K[iel, :] contains flattened element matrices for both DOF blocks:
            K[iel, :] = [Ael_A.flat, Ael_B.flat]
    sdof : int
        Total system degrees of freedom (2 * nnod for 2-DOF system)

    Returns
    -------
    A_all : scipy.sparse.csr_matrix, shape (sdof, sdof)
        Assembled global sparse matrix in CSR format.
        Ready for use with sparse linear solvers.

    Notes
    -----
    Assembly Process:
        1. Flatten all element contributions (I, J, K)
        2. Create COO sparse matrix (efficient for construction)
        3. Convert to CSR format (efficient for solving)

    Duplicate entries at the same (i, j) location are automatically summed,
    which is the correct behavior for FEM assembly.

    Examples
    --------
    >>> # After element assembly
    >>> I, J = create_assembly_indices_2dof(EL2NOD, nnodel=3)
    >>> # ... compute element matrices Ael_A, Ael_B ...
    >>> K[iel, :] = np.concatenate([Ael_A.flat, Ael_B.flat])
    >>> # Assemble global matrix
    >>> A_all = assemble_sparse_matrix_2dof(I, J, K, sdof=2*nnod)

    See Also
    --------
    create_assembly_indices_2dof : Generate assembly indices
    scipy.sparse.csr_matrix : CSR sparse matrix format
    scipy.sparse.linalg.spsolve : Solve sparse linear system
    """
    # Flatten arrays to 1D for COO format
    # Shape changes: (nel, entries_per_element) → (nel * entries_per_element,)
    I_flat = I.reshape(-1)
    J_flat = J.reshape(-1)
    K_flat = K.reshape(-1)

    # Create sparse matrix in COO format, then convert to CSR
    # COO format is efficient for construction
    # CSR format is efficient for linear algebra operations
    A_all = csr_matrix((K_flat, (I_flat, J_flat)), shape=(sdof, sdof))

    return A_all


def get_triangle_quadrature(nip):
    """
    Get Gauss quadrature points and weights for triangular elements.

    This is a convenience wrapper that provides integration rules for
    triangular elements. Currently supports 3-point quadrature (degree 2
    precision), which is sufficient for linear elements.

    Parameters
    ----------
    nip : int
        Number of integration points. Currently only nip=3 is supported.
        For more quadrature rules, use fem_integration.ip_triangle.

    Returns
    -------
    gauss : ndarray, shape (nip, 2)
        Integration point coordinates in reference element (xi, eta).
        Each row is one integration point.
    weights : ndarray, shape (nip,)
        Integration weights.
        sum(weights) = 0.5 (area of reference triangle)

    Notes
    -----
    For triangular elements with vertices at (0,0), (1,0), (0,1),
    the 3-point Gauss quadrature rule integrates polynomials up to
    degree 2 exactly.

    The integration formula is:
        ∫∫_T f(x,y) dA ≈ sum_i w_i * f(xi_i, eta_i) * detJ

    where w_i are weights, (xi_i, eta_i) are integration points,
    and detJ is the Jacobian determinant (= 2 * element_area).

    Examples
    --------
    >>> gauss, weights = get_triangle_quadrature(nip=3)
    >>> print("Integration points:")
    >>> print(gauss)
    [[0.16666667 0.16666667]
     [0.66666667 0.16666667]
     [0.16666667 0.66666667]]
    >>> print(f"Sum of weights: {np.sum(weights):.4f}")
    Sum of weights: 0.5000

    See Also
    --------
    fem_integration.ip_triangle : More comprehensive quadrature rules
    """
    if nip == 3:
        # 3-point Gauss quadrature for triangles (degree 2 precision)
        # Points located at 1/6 and 2/3 along edges
        gauss = np.array([
            [1/6, 1/6],   # Near vertex 1
            [2/3, 1/6],   # Near vertex 2
            [1/6, 2/3]    # Near vertex 3
        ])
        weights = np.array([1/6, 1/6, 1/6])
    else:
        raise NotImplementedError(
            f"Only nip=3 is currently supported in this function. "
            f"For nip={nip}, use fem_integration.ip_triangle instead."
        )

    return gauss, weights


if __name__ == "__main__":
    """
    Demonstration and testing of fem_utils functions.
    """
    print("=" * 70)
    print("FEM Utilities - Demonstration")
    print("=" * 70)
    print()

    # Test 1: Jacobian computation
    print("Test 1: Jacobian computation for unit right triangle")
    print("-" * 70)
    ECOORD = np.array([[0, 0], [1, 0], [0, 1]])
    dNds = np.array([[-1, 1, 0], [-1, 0, 1]])  # Linear element derivatives
    Jac, invJ, detJ = compute_jacobian(dNds, ECOORD)
    print(f"Jacobian matrix:\n{Jac}")
    print(f"Inverse Jacobian:\n{invJ}")
    print(f"Determinant: {detJ:.4f}")
    print(f"Element area: {detJ/2:.4f} (should be 0.5 for unit right triangle)")
    print()

    # Test 2: Assembly indices
    print("Test 2: Assembly indices for 2-DOF system")
    print("-" * 70)
    EL2NOD = np.array([[0, 1, 2], [1, 3, 2]])  # 2 elements
    I, J = create_assembly_indices_2dof(EL2NOD, nnodel=3)
    print(f"Element 0 connects nodes {EL2NOD[0]}")
    print(f"  A-DOF indices: {I[0, :9]}")
    print(f"  B-DOF indices: {I[0, 9:]}")
    print()

    # Test 3: Quadrature rule
    print("Test 3: Triangle quadrature rule")
    print("-" * 70)
    gauss, weights = get_triangle_quadrature(nip=3)
    print(f"Integration points:\n{gauss}")
    print(f"Weights: {weights}")
    print(f"Sum of weights: {np.sum(weights):.6f} (should be 0.5)")
    print()

    print("=" * 70)
    print("All demonstrations complete!")
    print("=" * 70)
