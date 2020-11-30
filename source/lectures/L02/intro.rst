Navier-Stokes flow with OpenFoam
================================

Theoretical background
-------------------------------------
The dynamics of fluids is described by the Navier-Stokes equation. It is based on mass and momentum conservation in a moving fluid.

For simplified incompressible case, these equations take the following form:

.. math::
    :label: eq:continuity 

    \nabla \cdot \vec{U} = 0
    
.. math::
    :label: eq:mom_con
    
    \frac{\partial \vec{U} }{\partial t} + \nabla \cdot (\vec{U} \vec{U}) = \nabla \cdot (\nu \nabla \vec{U}) - \nabla p

Solution strategies - the SIMPLE algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The 3-D Navier-Stokes equation has four unknowns (three velocity components and one pressure) but the standard formulation does not include an explicit "pressure equation", which complicated solving the system of equations.

.. tip::
    There is a very nice `youtube video  <https://www.youtube.com/watch?v=OOILoJ1zuiw>`_ by FluidDynamics101 (Aidan Winshurst) that explains the SIMPLE algorithm. We are here following his derivation and notation.

There are different solution strategies implemented in OpenFoam. The simplest one is the SIMPLE algorithm for the steady-state form of the incompressible Navier-Stokes equation.

.. math::
    :label: eq:continuity2 

    \nabla \cdot \vec{U} = 0
    
.. math::
    :label: eq:mom_con2
    
    \nabla \cdot (\vec{U} \vec{U}) = \nabla \cdot (\nu \nabla \vec{U}) - \nabla p

The SIMPLE algorithm solves these equations by deriving a pressure equation by combining the momentum and continuity equations and by formulating a corrector equation to make the velocity field satisfy the continuity equation.

Step 1 is to formulate the momentum balance as general matrix equation:


.. math::
    :label: eq:mom_con3
    
    \nabla \cdot (\vec{U} \vec{U}) = \nabla \cdot (\nu \nabla \vec{U}) - \nabla p

.. math::
    :label: eq:mom_matrix 

    M \vec{U} = - \nabla p
    
Next we split the matrix :math:`M` into a diagonal part and remaining part. 

.. math::
    :label: eq:mom_matrix2 

    M \vec{U} = - \nabla p

.. math::
    :label: eq:mom_matrix3

    A \vec{U} -H = - \nabla p

The nice thing about diagonal matrices is that they are easy to invert.

Step 2 is to solve for :math:`\vec U` using an initial guess for the pressure field.

.. math::
    :label: eq:mom_matrix4

    \vec{U} = A^{-1} H - A^{-1}A \nabla p

and plug it into the continuity equation:

.. math::
    :label: eq:mom_matrix5

        \nabla \cdot \vec{U} = 0
 
.. math::
    :label: eq:mom_matrix6      
    
        \nabla \cdot (A^{-1} H - A^{-1}A \nabla p) = 0

to get

.. math::
    :label: eq:mom_matrix7

        \nabla \cdot A^{-1} \nabla p = \nabla \cdot(A^{-1}H)

Equations :eq:`eq:mom_matrix` and :eq:`eq:mom_matrix6` form a set of equations that are solved in sequence.

First solve the momentum balance using an initial guess:

.. math::
    :label: eq:mom_matrix8 

    M \vec{U} = - \nabla p

Now update pressure by solving the pressure equation :eq:`eq:mom_matrix6`

.. math::
    :label: eq:mom_matrix9

        \nabla \cdot A^{-1} \nabla p = \nabla \cdot(A^{-1}H)

Use the pressure to "correct" the velocity solution, so that it fulfills the continuity equation :eq:`eq:continuity2`.

.. math::
    :label: eq:u_cor

        \vec{U} = A^{-1} H - A^{-1}A \nabla p

This is the SIMPLE algorithm for the steady-state incompressible Navier-Stokes equation. Later we will also learn about other algorithms for more general forms of the N-S stokes like the PISO and PIMPLE algorithms.