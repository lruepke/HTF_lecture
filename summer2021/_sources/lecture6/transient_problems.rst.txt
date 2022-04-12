Transient heat diffusion
=========================================
The next step is to derive and implement the unsteady (time-dependent) heat diffusion equation and to solve an example problem for the cooling of the lithosphere. The unsteady heat diffusion equation looks like this:

.. math::
    :label: eq:fem_2d_transient

    \rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} k \frac{\partial T}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T}{\partial y}
    
Two differences to the previously considered steady heat diffusion are apparent. First, we now have a time derivative in addition to the spatial derivatives. Second, the material parameters (density, specific heat, thermal conductivity) are not constant (neglected) any more. Nevertheless, the diffusion part of :eq:`eq:fem_2d_transient` looks very similar to the steady-state case that we know how to solve.

We know how to handle spatial derivatives but this is the first time we encounter a time derivative in finite elements. To solve it we will use a “trick” from finite differences – we will simply write the time derivative in finite differences form.

.. math::
    :label: eq:fem_2d_transient_2

    \rho c_p \frac{T^{n+1} - T^n}{\Delta t} = \frac{\partial}{\partial x} k\frac{\partial T^{n+1}}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T^{n+1}}{\partial y}.

Re-arrange :eq:`eq:fem_2d_transient` so that all known temperatures :math:`T^n` are on the Rhs and all unknown temperatures :math:`T^{n+1}` are on the Lhs.


.. math::
    :label: eq:fem_2d_transient_3

    \rho c_p T^{n+1} - \Delta t \left( \frac{\partial}{\partial x} k\frac{\partial T^{n+1}}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T^{n+1}}{\partial y} \right ) = \rho c_p T^{n}.


FEM form
--------

Now we proceed in the usual way: insert the approximate solution using shape functions and use the Galerkin method:

.. math::
    :label: eq:fem_2d_transient_weak

    \int_\Omega  \rho c_p N_i N_j T^{n+1}_j d\Omega -  \int_\Omega N_i  \Delta t \left ( \frac{\partial}{\partial x}k\frac{\partial N_j T^{n+1}_j }{\partial x} + N_i \frac{\partial}{\partial y}k\frac{\partial N_j T^{n+1}_j }{\partial y} \right ) d\Omega= \int_\Omega  \rho c_p N_i N_j T^{n}_j d\Omega   0\ \ \ \ \ \ \ i=1,2,...,n


and integrate the diffusion term by parts:

.. math::
    :label: eq:fem_2d_transient_weak_v2

    \int_\Omega  \rho c_p N_i N_j T^{n+1}_j d\Omega +  \int_\Omega \left ( \Delta t  \frac{\partial N_i}{\partial x}k\frac{\partial N_j T_j }{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j T_j }{\partial y} d\Omega \right ) d\Omega= \int_\Omega  \rho c_p N_i N_j T^{n}_j d\Omega  - \oint_{\Gamma} N_i \Delta t \vec{q}\vec{n}  d\Gamma\ \ \ \ \ \ \ i=1,2,...,n


We can proceed and write everything in terms of matrices:

.. math::
    :label: eq:fem_2d_transient_weak_matrix

    \left( \rho c_p M  + \Delta t A  \right ) T^{n+1} =  \rho c_p M  T^{n} + BC \\

With the matrices being define as:

.. math::
    :label: eq:fem_2d_transient_weak_matrix_v2

    \begin{align}
    \begin{split}
    M &= \int_\Omega N_i N_j  d\Omega \ \ \ \ \ \ \ i,j=1,2,...,n\\
    A &= \int_\Omega \left ( \frac{\partial N_i}{\partial x}k\frac{\partial N_j }{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j }{\partial y} d\Omega \right ) d\Omega\ \ \ \ \ \ \ i,j=1,2,...,n\\
    \end{split}
    \end{align}


The matrix :math:`M` is called the mass matrix. The terms in brackets on the LHS of :eq:`eq:fem_2d_transient_weak_matrix` will become the new matrix that is assembled per element and added to the global stiffness matrix.

Implementation
--------------

We implement the transient behavior into our triangle script from the previous lecture. If you didn't complete it, you can download it from here (:download:`2d_fem_transient_triangle.py <python/2d_fem_transient_triangle.py>`).

We will have to make several changes to the code:

Time loop and output writing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For transient problems, we will need a time loop over all time steps and a way to store/visualize the evolving solution. So far, we have plotted only steady-state solutions that we could directly plot using :code:`matplotlib`. Now we will have to come up with a different strategy as we want to visualize and analyze the complete transient solution. One good way is to use the :code:`meshio` python package to save the solution in vtk format and to later analyze is using `Paraview <https://www.paraview.org>`_ If you don't have paraview installed, no is the time ;)

Here is a code snippit that illustrates the logic:

.. code-block:: python 
    :caption: Output writing.

    # mesh information
    points=np.hstack((GCOORD, GCOORD[:,0].reshape(-1,1)*0)) #must have 3 components (x,y,z)
    cells=[("triangle",EL2NOD)]
    
    # write initial mesh
    writer=meshio.xdmf.TimeSeriesWriter('transient.xmf')
    writer.__enter__() # have to add this: import hdf5 and open file ...
    writer.write_points_cells(points, cells)
    
    dt = 0.025
    nt = 40

    # model time loop
    for t in range(0,nt):
        #our FEM code here

        #save results
        #cell data
        U=np.hstack((Q_x.reshape(-1,1),Q_y.reshape(-1,1)))
        U=np.hstack((U,U[:,0].reshape(-1,1)*0))
    
        #save data
        writer.write_data(t, point_data={"T": T},cell_data={"U": [U], "K": [Kel]})
    writer.__exit__() # close file

Note the time loop over all time steps. The number of time steps and the time step itself are chosen to make the final results look nice - they don't have a physical meaning for the time being.

Note also, that we don't need the python plotting at the end of the script anymore.


Matrix assembly
^^^^^^^^^^^^^^^
If we look at :eq:`eq:fem_2d_transient_weak_matrix_v2`, we notice that we have to change the matrix assembly to 1) account for the mass matrix in the element stiffness matrix, and 2) to integrate the old temperatures into the force vector. This can be done like this:


.. code-block:: python 
    :caption: matrix assembly.

    # 4. compute element stiffness matrix
    Ael     = Ael + (rho*cp*np.outer(N,N) +  dt*Kel[iel]*np.matmul(dNdx.T, dNdx))*detJ*weights[ip] # [nnodel,1]*[1,nnodel] / weights are missing, they are 1
    
    # 5. assemble right-hand side
    Rhs_el     = Rhs_el + rho*cp*np.matmul(np.outer(N,N), np.take(T, EL2NOD[iel,:], axis=0 ))*detJ*weights[ip] 


Notice how the logic for the element thermal conductivity has changed - and that we need two additional physical parameters :math:`\rho` and :math:`c_p` . 

.. code-block:: python 
    :caption: model parameters.

    rho         = 1
    cp          = 1

    Kel    = np.ones(nel)*k1
    Kel[np.where(Phases==100)] = k2

Make sure that the new logical for :code:`Kel` is also used in the post-processing step when computing heat fluxes.

.. only:: html

    Results of transient diffusion problem.
 
    .. raw:: html
 
       <video width=100% autoplay muted controls loop>
       <source src="../_static/video/T.mp4" type="video/mp4">
          Your browser does not support HTML video.
       </video>


