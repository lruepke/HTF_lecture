.. _blocked_assembly:

Blocked Vectorized Assembly
===========================

In the previous sections, we implemented a 2D Stokes flow solver using an element-by-element assembly loop. While correct and readable, this approach becomes a performance bottleneck for large meshes. In this section, we explain how to refactor the assembly loop using a **blocked vectorized** strategy that processes many elements simultaneously using NumPy array operations.

The key idea is simple: instead of looping over elements one at a time, we group elements into **blocks** and perform all operations on the entire block at once using NumPy's vectorized operations and ``np.einsum``.


Why element-by-element assembly is slow
-----------------------------------------

In our original implementation (``mechanical2d.py``), the assembly loop processes one element at a time:

.. code-block:: python

    for iel in range(nel):
        ECOORD_X = mesh.GCOORD[mesh.EL2NOD[iel, :], :]
        EMu = material.Mu[mesh.Phases[iel]]

        A_elem = np.zeros((nedof, nedof))
        # ... integration point loop ...
        for ip in range(solver.nip):
            Jac = dNdui @ ECOORD_X
            invJ = np.linalg.inv(Jac)
            # ... build B, accumulate A_elem, Q_elem, M_elem ...

        # Static condensation
        invM_elem = np.linalg.inv(M_elem)
        A_elem += PF * (Q_elem.T @ invM_elem @ Q_elem)

        A_all[iel, :] = A_elem.ravel()

This has two problems:

1. **Python loop overhead**: Each iteration through a Python ``for`` loop incurs interpreter overhead. For a mesh with 10,000 elements and 6 integration points, that is 60,000 iterations of pure Python.

2. **Poor cache utilization**: Processing one element at a time means each ``np.linalg.inv`` or matrix multiply operates on a single tiny matrix (2x2 or 3x3). The CPU spends more time on function call overhead than on actual computation.

NumPy is fast when it operates on **large arrays**. A single call to invert 10,000 matrices at once is much faster than 10,000 calls to invert one matrix each.


The blocking strategy
-----------------------

The blocked approach divides elements into chunks (blocks) and processes each block using vectorized NumPy operations:

.. code-block:: text

    Elements:  [0, 1, 2, ..., 9999, 10000, 10001, ..., 19999, 20000, ...]
               |---- Block 0 -------|------- Block 1 ---------|-- Block 2 ...

    For each block:
        - Extract coordinates for ALL elements in block     (n_el, 7, 2)
        - For each integration point (small loop, 6 iterations):
            - Compute ALL Jacobians at once                  (n_el, 2, 2)
            - Invert ALL Jacobians at once                   (n_el, 2, 2)
            - Build ALL B matrices at once                   (n_el, 3, 14)
            - Accumulate ALL element matrices at once        (n_el, 14, 14)
        - Static condensation for ALL elements at once

The block size is a tuning parameter. A block of 10,000 elements uses roughly:

.. math::

    10{,}000 \times 14 \times 14 \times 8\,\text{bytes} \approx 15\,\text{MB}

This fits comfortably in modern CPU L3 caches (typically 8--32 MB), which means data stays close to the CPU during the integration point loop.


Extracting block data
-----------------------

The first step in each block is to gather coordinates and material properties for all elements at once:

.. code-block:: python

    for block_idx in range(n_blocks):
        start_el = block_idx * block_size
        end_el = min((block_idx + 1) * block_size, nel)
        n_el_block = end_el - start_el

        # Fancy indexing: extract coordinates for all elements in block
        EL2NOD_block = mesh.EL2NOD[start_el:end_el]         # (n_el, 7)
        ECOORD_block = mesh.GCOORD[EL2NOD_block]            # (n_el, 7, 2)
        EMu_block = material.Mu[mesh.Phases[start_el:end_el]]  # (n_el,)

The expression ``mesh.GCOORD[EL2NOD_block]`` uses NumPy's advanced indexing: ``EL2NOD_block`` has shape ``(n_el, 7)`` containing node indices, and indexing into ``GCOORD`` (shape ``(nnod, 2)``) returns an array of shape ``(n_el, 7, 2)`` --- the x,y coordinates of all 7 nodes for every element in the block.


Vectorized Jacobian computation
---------------------------------

In the original code, the Jacobian is computed for one element:

.. code-block:: python

    # Original: single element, shapes (2, 7) @ (7, 2) = (2, 2)
    Jac = dNdui @ ECOORD_X
    invJ = np.linalg.inv(Jac)

In the blocked version, we compute Jacobians for the entire block using ``np.einsum``:

.. code-block:: python

    # Blocked: all elements at once
    # dNdui:        (2, 7)       - same for all elements (reference element)
    # ECOORD_block: (n_el, 7, 2) - different per element
    # Result:       (n_el, 2, 2) - one Jacobian per element
    Jac_block = np.einsum('ik,nkj->nij', dNdui, ECOORD_block)


Understanding ``np.einsum``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Einstein summation convention ``'ik,nkj->nij'`` reads as:

- ``ik`` --- the first array has indices ``i`` (rows) and ``k`` (columns)
- ``nkj`` --- the second array has indices ``n`` (element), ``k`` (node), ``j`` (x/y)
- ``nij`` --- the output has indices ``n`` (element), ``i`` (row), ``j`` (col)
- The repeated index ``k`` is summed over (matrix contraction)

This is equivalent to ``Jac[n, i, j] = sum_k dNdui[i, k] * ECOORD[n, k, j]`` for every element ``n`` --- a batched matrix multiply.


Vectorized 2x2 matrix inversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rather than calling ``np.linalg.inv`` on each 2x2 Jacobian, we invert all of them using the analytical formula:

.. code-block:: python

    def batch_invert_2x2(matrices):
        """Invert n matrices of shape (n, 2, 2) simultaneously."""
        a = matrices[:, 0, 0]
        b = matrices[:, 0, 1]
        c = matrices[:, 1, 0]
        d = matrices[:, 1, 1]

        det = a * d - b * c

        inv = np.zeros_like(matrices)
        inv[:, 0, 0] =  d / det
        inv[:, 0, 1] = -b / det
        inv[:, 1, 0] = -c / det
        inv[:, 1, 1] =  a / det

        return inv, det

Each line operates on arrays of length ``n_el``, so inverting 10,000 Jacobians is essentially as fast as inverting one.


Vectorized B matrix and stiffness assembly
---------------------------------------------

The strain-displacement matrix :math:`\mathbf{B}` maps nodal velocities to strain rates. In the original code, it is built for one element:

.. code-block:: python

    # Original: single element
    B[0, 0::2] = dNdX[0]  # dN/dx for x-velocity DOFs
    B[1, 1::2] = dNdX[1]  # dN/dy for y-velocity DOFs
    B[2, 0::2] = dNdX[1]  # dN/dy for x-velocity DOFs (shear)
    B[2, 1::2] = dNdX[0]  # dN/dx for y-velocity DOFs (shear)

In the blocked version, this becomes:

.. code-block:: python

    # Blocked: all elements, B_block has shape (n_el, 3, 14)
    B_block = np.zeros((n_el_block, 3, nedof))
    B_block[:, 0, 0::2] = dNdX_block[:, 0, :]   # (n_el, 7)
    B_block[:, 1, 1::2] = dNdX_block[:, 1, :]
    B_block[:, 2, 0::2] = dNdX_block[:, 1, :]
    B_block[:, 2, 1::2] = dNdX_block[:, 0, :]

The element stiffness :math:`\mathbf{A}^e = \int \eta\, \mathbf{B}^T \mathbf{D}\, \mathbf{B}\, d\Omega` is then computed for all elements at once:

.. code-block:: python

    # DEV is (3, 3) - constant for all elements
    # B_block is (n_el, 3, 14)

    # Step 1: D @ B for all elements → (n_el, 3, 14)
    DEV_B = np.einsum('ij,njk->nik', DEV, B_block)

    # Step 2: B^T @ (D @ B) for all elements → (n_el, 14, 14)
    BtDB_block = np.einsum('nji,njk->nik', B_block, DEV_B)

    # Step 3: Scale by weight * viscosity per element
    A_block += (weight_block * EMu_block)[:, np.newaxis, np.newaxis] * BtDB_block

The ``[:, np.newaxis, np.newaxis]`` broadcasting trick turns the per-element scalar ``weight * mu`` from shape ``(n_el,)`` into shape ``(n_el, 1, 1)`` so it can multiply the ``(n_el, 14, 14)`` stiffness matrices element-wise.


Vectorized pressure shape functions
--------------------------------------

The pressure interpolation matrix :math:`\mathbf{P}` and vector :math:`\mathbf{P}_b` are element-specific because they depend on the corner node coordinates. In the original code:

.. code-block:: python

    # Original: per element
    P = np.ones((3, 3))
    P[1:3, :] = ECOORD_X[:3].T     # corner coordinates
    # Inside ip loop:
    Pb[1:3] = Ni @ ECOORD_X        # integration point coordinates
    Pi = np.linalg.solve(P, Pb)    # pressure shape functions

In the blocked version, we construct :math:`\mathbf{P}` for all elements and invert it once per block (outside the integration point loop):

.. code-block:: python

    # Blocked: P for all elements, shape (n_el, 3, 3)
    P_block = np.ones((n_el_block, 3, 3))
    P_block[:, 1:3, :] = ECOORD_block[:, :3, :].transpose(0, 2, 1)

    # Invert P once per block (reused for all integration points)
    P_inv_block = np.linalg.inv(P_block)   # (n_el, 3, 3)

    # Inside the ip loop: compute Pi = P^{-1} @ Pb for all elements
    ip_coords = np.einsum('k,nkj->nj', Ni, ECOORD_block)  # (n_el, 2)
    Pb_block = np.ones((n_el_block, 3))
    Pb_block[:, 1] = ip_coords[:, 0]
    Pb_block[:, 2] = ip_coords[:, 1]

    Pi_block = np.einsum('nij,nj->ni', P_inv_block, Pb_block)  # (n_el, 3)

Since :math:`\mathbf{P}` only depends on corner coordinates (constant per element), precomputing :math:`\mathbf{P}^{-1}` avoids calling ``np.linalg.solve`` at every integration point. The batched ``np.linalg.inv`` handles all elements in one call.


Vectorized pressure coupling and mass matrix
-----------------------------------------------

The pressure-velocity coupling :math:`\mathbf{Q}` and pressure mass matrix :math:`\mathbf{M}` involve outer products of the pressure shape functions :math:`\Pi` with the volumetric strain operator :math:`\mathbf{B}_{\text{vol}}`:

.. code-block:: python

    # Original: per element
    Q_elem -= weight * np.outer(Pi, Bvol.ravel())     # (3,) x (14,) → (3, 14)
    M_elem += weight * np.outer(Pi, Pi)                # (3,) x (3,)  → (3, 3)

In the blocked version, ``np.einsum`` computes outer products for all elements simultaneously:

.. code-block:: python

    # Blocked: outer products for all elements
    # Pi_block: (n_el, 3), Bvol_ravel_block: (n_el, 14)
    Q_block -= weight_block[:, None, None] * np.einsum('ni,nj->nij', Pi_block, Bvol_ravel_block)
    M_block += weight_block[:, None, None] * np.einsum('ni,nj->nij', Pi_block, Pi_block)

The volumetric strain :math:`\mathbf{B}_{\text{vol}}` is obtained by rearranging the global shape function derivatives:

.. code-block:: python

    # dNdX_block has shape (n_el, 2, 7)
    # We need [dNdx_0, dNdy_0, dNdx_1, dNdy_1, ...] per element → (n_el, 14)
    Bvol_ravel_block = dNdX_block.transpose(0, 2, 1).reshape(n_el_block, -1)


Vectorized static condensation
---------------------------------

After the integration point loop, the pressure DOFs are eliminated via static condensation. The original code does this per element:

.. code-block:: python

    # Original: per element
    invM_elem = np.linalg.inv(M_elem)                      # (3, 3)
    A_elem += PF * (Q_elem.T @ invM_elem @ Q_elem)         # (14, 14)

The blocked version uses batched inversion and ``np.einsum`` for the triple matrix product:

.. code-block:: python

    # Blocked: static condensation for all elements
    invM_block = np.linalg.inv(M_block)                         # (n_el, 3, 3)
    invM_Q = np.einsum('nij,njk->nik', invM_block, Q_block)    # (n_el, 3, 14)
    QtinvMQ = np.einsum('nji,njk->nik', Q_block, invM_Q)       # (n_el, 14, 14)
    A_block += PF * QtinvMQ

Here ``np.linalg.inv`` natively supports batch dimensions, so inverting all 3x3 pressure mass matrices is a single function call.


Scattering the right-hand side
---------------------------------

The right-hand side vector requires scattering element contributions to global DOFs. Since multiple elements share nodes, we cannot simply use array indexing (which would silently drop duplicate contributions). Instead, ``np.add.at`` performs unbuffered in-place addition:

.. code-block:: python

    # Scatter RHS from block to global vector
    for i in range(nedof):
        np.add.at(Rhs_all, EL2DOF_block[:, i], Rhs_block[:, i])

This replaces the original ``Rhs_all[EL2DOF[iel, :]] += Rhs_elem`` which only handled one element at a time. The small loop over ``nedof`` (14 iterations) is unavoidable because ``np.add.at`` needs a 1D index array.


Sparse matrix index construction
-----------------------------------

In the original code, sparse matrix row/column indices are built inside the element loop:

.. code-block:: python

    # Original: per element
    I[iel, :] = np.tile(EL2DOF[iel, :], (nedof, 1)).T.ravel()
    J[iel, :] = np.tile(EL2DOF[iel, :], (nedof, 1)).ravel()

In the blocked version, this is done once for **all** elements outside the block loop:

.. code-block:: python

    # Vectorized: all elements at once
    I_all = np.repeat(EL2DOF, nedof, axis=1)   # each DOF repeated nedof times
    J_all = np.tile(EL2DOF, (1, nedof))         # DOF array tiled nedof times

``np.repeat`` and ``np.tile`` operate on the full ``(nel, 14)`` array, producing the ``(nel, 196)`` index arrays in one call.


Summary of array shapes
--------------------------

The following table summarizes how array dimensions change from element-by-element to blocked processing:

.. list-table:: Array shapes: original vs. blocked
   :widths: 30 25 25
   :header-rows: 1

   * - Quantity
     - Original (per element)
     - Blocked (per block)
   * - Coordinates
     - ``(7, 2)``
     - ``(n_el, 7, 2)``
   * - Jacobian
     - ``(2, 2)``
     - ``(n_el, 2, 2)``
   * - Shape function derivatives
     - ``(2, 7)``
     - ``(n_el, 2, 7)``
   * - B matrix
     - ``(3, 14)``
     - ``(n_el, 3, 14)``
   * - Element stiffness A
     - ``(14, 14)``
     - ``(n_el, 14, 14)``
   * - Pressure coupling Q
     - ``(3, 14)``
     - ``(n_el, 3, 14)``
   * - Pressure mass M
     - ``(3, 3)``
     - ``(n_el, 3, 3)``
   * - Pressure shape functions Pi
     - ``(3,)``
     - ``(n_el, 3)``

The pattern is consistent: every per-element quantity gains an extra leading dimension ``n_el`` (the block size).


Performance comparison
-------------------------

For a mesh with approximately 9,500 elements (7-node triangles, 6 integration points):

.. list-table:: Assembly time comparison
   :widths: 35 20 20
   :header-rows: 1

   * - Version
     - Assembly time
     - Speedup
   * - Original (element-by-element)
     - ~100 ms
     - 1.0x
   * - Blocked vectorized
     - ~5 ms
     - ~20x

The overall solver speedup is smaller because the linear solve (``spsolve``) dominates for large systems, but the assembly phase --- previously the bottleneck for moderately-sized meshes --- is effectively eliminated as a concern.

The blocked implementation has been verified to produce results identical to the original (velocity differences < :math:`10^{-12}`, pressure differences < :math:`10^{-10}`).


Download
---------

You can download the blocked implementation here: :download:`mechanical2d_blocked.py <python/mechanical2d_blocked.py>` and :download:`mechanical2d_driver_blocked.py <python/mechanical2d_driver_blocked.py>`.

Exercises
---------

1. **Block size tuning**: Run the blocked solver with different ``BLOCK_SIZE`` values (100, 1000, 10000, 50000). How does assembly time change? Can you explain the behavior in terms of CPU cache sizes?

2. **Profiling**: Compare the timing breakdown between the original and blocked versions. At what mesh size does the assembly cease to be the bottleneck?

3. **einsum practice**: Write out the index notation for each ``np.einsum`` call in the blocked assembly. Verify on paper that the contraction indices are correct for one of the operations (e.g., the ``B^T D B`` computation).
