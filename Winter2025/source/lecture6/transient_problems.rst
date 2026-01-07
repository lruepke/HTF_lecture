Transient heat diffusion
=========================================
The next step is to derive and implement the unsteady (time-dependent) heat diffusion equation and to solve an example problem for the cooling of the lithosphere. The transient heat equation looks like this:

.. math::
    :label: eq:fem_2d_transient

    \rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} k \frac{\partial T}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T}{\partial y}
    
Two differences to the previously considered steady-state heat diffusion are apparent. First, we now have a time derivative in addition to the spatial derivatives. Second, the material parameters (density, specific heat, thermal conductivity) are no longer assumed constant. Nevertheless, the diffusion part of :eq:`eq:fem_2d_transient` looks very similar to the steady-state case that we know how to solve.

We know how to handle spatial derivatives but this is the first time we encounter a time derivative in finite elements. We use a first-order implicit (backward Euler) finite-difference for time:

.. math::
    :label: eq:fem_2d_transient_2

    \rho c_p \frac{T^{n+1} - T^n}{\Delta t} = \frac{\partial}{\partial x} k\frac{\partial T^{n+1}}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T^{n+1}}{\partial y}.

Re-arrange :eq:`eq:fem_2d_transient_2` so that all known temperatures :math:`T^n` are on the Rhs and all unknown temperatures :math:`T^{n+1}` are on the Lhs.


.. math::
    :label: eq:fem_2d_transient_3

    \rho c_p T^{n+1} - \Delta t \left( \frac{\partial}{\partial x} k\frac{\partial T^{n+1}}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T^{n+1}}{\partial y} \right ) = \rho c_p T^{n}.

The backward Euler scheme is unconditionally stable for linear diffusion (useful for larger time steps), but only first-order accurate in time.


FEM form
--------

Now we proceed in the usual way: insert the approximate solution using shape functions and use the Galerkin method. Writing :math:`T \approx \sum_j N_j T_j` and testing with :math:`N_i` gives the weak form:

.. math::
    :label: eq:fem_2d_transient_weak

    \int_\Omega  \rho c_p\, N_i N_j\, T^{n+1}_j\, d\Omega 
    +  \Delta t \int_\Omega \nabla N_i \cdot \big( k\, \nabla N_j \big)\, T^{n+1}_j\, d\Omega
    = \int_\Omega  \rho c_p\, N_i N_j\, T^{n}_j\, d\Omega  
    - \Delta t \oint_{\Gamma_N} N_i\, \vec{q}\cdot\vec{n}\, d\Gamma\ \ \ \ \ \ \ i=1,2,...,n


and we obtain this by integrating the diffusion term by parts (moving derivatives from :math:`T` to the test function) and collecting the Neumann boundary term on :math:`\Gamma_N`.

.. math::
    :label: eq:fem_2d_transient_weak_v2

    \int_\Omega  \rho c_p\, N_i N_j\, T^{n+1}_j\, d\Omega 
    +  \Delta t \int_\Omega \nabla N_i \cdot \big( k\, \nabla N_j \big)\, T^{n+1}_j\, d\Omega
    = \int_\Omega  \rho c_p\, N_i N_j\, T^{n}_j\, d\Omega  
    - \Delta t \oint_{\Gamma_N} N_i\, \vec{q}\cdot\vec{n}\, d\Gamma\ \ \ \ \ \ \ i=1,2,...,n


We can proceed and write everything in terms of matrices:

.. math::
    :label: eq:fem_2d_transient_weak_matrix

    \left( \rho c_p M  + \Delta t A  \right ) T^{n+1} =  \rho c_p M  T^{n} + BC \\

With the matrices defined as:

.. math::
    :label: eq:fem_2d_transient_weak_matrix_v2

    \begin{align}
    \begin{split}
    M &= \int_\Omega N_i N_j\, d\Omega \ \ \ \ \ \ \ i,j=1,2,...,n\\
    A &= \int_\Omega \left ( \frac{\partial N_i}{\partial x}\,k\,\frac{\partial N_j }{\partial x} + \frac{\partial N_i}{\partial y}\,k\,\frac{\partial N_j }{\partial y} \right ) d\Omega\ \ \ \ \ \ \ i,j=1,2,...,n\\
    \end{split}
    \end{align}


The matrix :math:`M` is called the mass matrix. The terms in brackets on the LHS of :eq:`eq:fem_2d_transient_weak_matrix` will become the new matrix that is assembled per element and added to the global stiffness matrix. If :math:`\rho c_p` or :math:`k` vary spatially, keep them inside the integrals (pulling them out assumes element-wise constants).

Implementation
--------------

We implement the transient behavior into our triangle script from the previous lecture. If you didn't complete it, you can download it from here (:download:`2d_fem_transient_triangle.py <python/2d_fem_transient_triangle.py>`).

We will have to make several changes to the code:

Time loop and output writing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For transient problems, we will need a time loop over all time steps and a way to store/visualize the evolving solution. So far, we have plotted only steady-state solutions that we could directly plot using :code:`matplotlib`. Now we will have to come up with a different strategy as we want to visualize and analyze the complete transient solution. One good way is to use the :code:`meshio` python package to save the solution in XDMF/HDF5 (VTK) format and then analyze it using `ParaView <https://www.paraview.org>`_. If you don't have ParaView installed, now is the time ;)

Here is a minimal code snippet that illustrates the logic:

.. code-block:: python 
    :caption: Output writing.

    import meshio

    # mesh information (x,y,0)
    points = np.hstack((GCOORD, np.zeros((GCOORD.shape[0], 1))))
    cells = [("triangle", EL2NOD)]

    dt = 0.025
    nt = 40

    with meshio.xdmf.TimeSeriesWriter('transient.xmf') as writer:
        writer.write_points_cells(points, cells)

        # model time loop
        for it in range(nt):
            # ... FEM update here ...

            # example cell data
            U = np.c_[Q_x, Q_y, np.zeros_like(Q_x)]
            writer.write_data(it*dt, point_data={"T": T}, cell_data={"U": [U], "K": [Kel]})

Note the time loop over all time steps. The number of time steps and the time step itself are chosen to make the final results look nice—they don't have a physical meaning for the time being.

Note also, that we don't need the python plotting at the end of the script anymore.


Matrix assembly
^^^^^^^^^^^^^^^
If we look at :eq:`eq:fem_2d_transient_weak_matrix_v2`, we notice that we have to change the matrix assembly to 1) account for the mass matrix in the element stiffness matrix, and 2) to integrate the old temperatures into the force vector. This can be done like this:


.. code-block:: python 
    :caption: matrix assembly.

    # 4. compute element matrix (mass + diffusion)
    Ael += (rho*cp*np.outer(N, N) + dt*Kel[iel]*(dNdx.T @ dNdx)) * detJ * weights[ip]

    # 5. assemble right-hand side from previous time step
    Rhs_el += rho*cp * (np.outer(N, N) @ T[EL2NOD[iel, :]]) * detJ * weights[ip]


Notice how the logic for the element thermal conductivity has changed - and that we need two additional physical parameters :math:`\rho` and :math:`c_p` . 

.. code-block:: python 
    :caption: model parameters.

    rho         = 1
    cp          = 1

    Kel    = np.ones(nel)*k1
    Kel[np.where(Phases==100)] = k2

Make sure that the new logic for :code:`Kel` is also used in the post-processing step when computing heat fluxes.

Remark (mass matrix): The expressions above use the consistent mass :math:`M_{ij}=\int N_i N_j\,d\Omega`. A lumped (diagonal) mass is also common and can simplify time stepping.

Remark (boundary conditions): Neumann heat fluxes contribute to the right-hand side via the boundary integral; Dirichlet temperatures are imposed as usual by modifying rows/columns.

.. only:: html

    Results of transient diffusion problem.
 
    .. raw:: html
 
       <video width=100% autoplay muted controls loop>
       <source src="../_static/video/T.mp4" type="video/mp4">
          Your browser does not support HTML video.
       </video>


