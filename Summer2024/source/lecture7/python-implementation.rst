.. _python-implementation:

Python implementation
===========================

We will use our previously developed 2d FEM code for steady-state diffusion as a starting point and modify it to solve for Stokes flow. There are a number of extra complications in the Stokes flow problem that we need to address:


Global degrees of freedom
------------------------------------------------------

In the steady-state diffusion problem, we had one degree of freedom per node. In the Stokes flow problem, we have two degrees of freedom per node: one for the x-velocity and one for the y-velocity. We will therefor have to double the number of degrees of freedom in the global stiffness matrix and the global force vector. We use a numberung where the 0th node has the degree of freedom 0 (x-velocity) and the degree of freedom 1 (y-velocity), the 1st node has the degrees of freedom 2 and 3, and so on. We will have to modify the code to account for this new numbering.

.. code-block:: python

    # MAKE ELDOF
    EL2DOF = np.zeros((nel,nedof), dtype=np.int32)
    EL2DOF[:,0::ndim] = ndim * EL2NOD
    EL2DOF[:,1::ndim] = ndim * EL2NOD + 1

We will use this matrix :math:`EL2DOF` to construct the global stiffness matrix and the global force vector.


Pressure shape functions
-------------------------

Step one is to spell out what the pressure shape functions are. In our derivation, we just stated that pressure is discontinuous and varies linearly over the element. We therefor have three unknown pressures per element. But what are those shape functions?

We start with the assumption that pressure varies linearly as a function of global x,y coordinates.

.. math::

    p(x,y) = a + bx + cy


Where :math:`a`,:math:`b`,:math:`c` are coefficients to be determined. The simplest linear shape function for a triangle with vertices at :math:`(x_1,y_1)`, :math:`(x_2,y_2)`, :math:`(x_3,y_3)` can be defined using barycentric coordinates or directly in terms of x and y coordinates.

We will find those shape functions by solving a linear system of equations. Let's define a matrix :math:`P` and a vector :math:`Pb` as follows:

.. math::

    P = \begin{bmatrix}
    1 & x_1 & y_1 \\
    1 & x_2 & y_2 \\
    1 & x_3 & y_3
    \end{bmatrix}

and a vector :math:`Pb` as follows:

.. math::
    Pb = \begin{bmatrix}
    1 \\
    x \\
    y
    \end{bmatrix}

where :math:`x` and :math:`y` are the coordinates of an integration point (or any point in the element). The interpolation functions (or shape functions) can now be found for solving the following linear system of equations:

.. math::

    P \begin{bmatrix}
    \Pi_1 \\
    \Pi_2 \\
    \Pi_3
    \end{bmatrix} = Pb

And we can find the shape functions as follows:

.. math::

   \Pi = P^{-1} Pb


In the python FEM implementation, we will spell out the matrix :math:`P` and the vector :math:`Pb` using the coordinates of the vertices of the element and the coordinates of the integration point; solving the matrix equation above will give us the pressure shape functions.

.. code-block:: python

        # np_edof is the number of pressure degrees of freedom (here 3)
        # ECOORD_X is a 7x2 matrix constructed ECOORD_X = GCOORD[EL2NOD[iel,:],:], holding the x and y coordinates of the element vertices
        P = np.ones((np_edof, np_edof))
        P[1:3, :] = ECOORD_X[:3].T


Inside the integratin loop, we construct :math:`Pb` and solve for the pressure shape functions:

.. code-block:: python

    for ip in range(nip):
        # Ni is [7,] and holds the velocity shape functions
        # Ni @ ECOORD_X is then [2,] and holds the x and y coordinates of the integration point
        # Pi is [3,] and holds the pressure shape functions

        Pb[1:3] = Ni @ ECOORD_X
        Pi      = np.linalg.solve(P, Pb)



Assembly of global matrices and pressure numbering
---------------------------------------------------

Finally, we need to assemble the global matrices. This includes the global :math:`Q` and :math:`invM` matrices, which we will use to reconstruct pressure from the velocity solution :ref:`eq:pressure_static_condensation` . 

As the pressure is discontinuous, each element has three independent pressure unknowns. The numbering is therefore that element 0 has unknown pressures 0,1,2 and element 1 has unknown pressures 3,4,5, and so on. For the global :math:`Q` matrix, the rows (i index) are therefore given by this continuous pressure numbering and the columns (j-index) by the matrix :math:`EL2DOF`.

.. code-block:: python

    # logic is that each element as three pressure dofs, as the pressure is dicsontinuous, we make ad-hoc numbering. 
    # element[0] -> pressure = [0,1,2]
    # element[1] -> pressure = [3,4,5]
    # The velocity dofs come from EL2DOF

    Q_i     = np.tile(np.arange(0, nel*np_edof, dtype=np.int32), (nedof,1)).T
    Q_j     = np.tile(EL2DOF, (1,np_edof))
    Q_all   = csr_matrix((Q_all.ravel(), (Q_i.ravel(), Q_j.ravel())), shape=(nel*np_edof, sdof))


The assembly of the global :math:`invM` matrix is similar. The rows and columns are both given by the continuous pressure numbering.

.. code-block:: python

    invM_i          = np.tile(np.arange(0, nel*np_edof, dtype=np.int32), (np_edof, 1)).T
    base_sequence   = np.tile(np.arange(np_edof), nel * np_edof)
    offsets         = np.repeat(np.arange(nel) * np_edof, np_edof**2)
    column_indices  = base_sequence + offsets
    invM_all        = csr_matrix((invM_all.ravel(), (invM_i.ravel(), column_indices.ravel())), shape=(nel*np_edof, nel*np_edof))    


There might be a more elegant solution to this; we'd be happy to hear about it!


Test problem
------------
We will use a test problem from structural geology to test our implementation. The problem is a pure shear problem with velocity boundary conditions at the sides and one (or multiple) inclusions of variable viscosity inside the modeling domain. Depending on the viscosity contrast, these inclusions will result in pressure anomalies. The problem has been explored analytically in :cite:`Schmid2003`. 

.. figure:: /_figures/pressure_field.png
    :name: fig:pressure_inclusion
    :align: center


You can download the code for this lecture here: :download:`ip_triangle.py <python/ip_triangle.py>`, :download:`shp_deriv_triangle.py <python/shp_deriv_triangle.py>`, :download:`mechanical2d_driver.py <python/mechanical2d_driver.py>`, and :download:`mechanical2d.py <python/mechanical2d.py>`.

Excerices
-----------
Get the code to work and try to solve the following exercises:

1. Modify the code to solve for a different problem, e.g., the pure shear problem with inclusions of different viscosity.
2. Modify the code to resolve inclusions of different shape
3. Modify the code to resolve inclusions of different densiy. This involves changes the boundary conditions to no slip. If you want, add a pseudo time loop in which GCOORD is updated using the computed velocities. This will give you a time-dependent problem with a moving inclusion.

