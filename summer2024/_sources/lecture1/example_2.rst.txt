Example 2: Cooling dike - implicit version
==========================================

In the previous session, we have learned how to solve for temperature using finite differences. However, our solution was only stable for certain parameter choices. In fact, we found that the finite differences formulation is only stable for :math:`\beta<0.5`:

.. math::
    :label: eq:imp_dummy1

    \begin{align}
    \begin{split}
    T_i^{n+1} = \frac{\kappa \Delta t}{\Delta x^2} \left( T_{i-1}^n - 2T_i^n + T_{i+1}^n\right) + T_i^n\\
    \beta = \frac{\kappa \Delta t}{\Delta x^2}\\
    T_i^{n+1} =
    \begin{bmatrix}
    \beta & (1-2\beta) & \beta
    \end{bmatrix}
    \begin{bmatrix}
    T_{i-1}^n\\
    T_{i}^n\\
    T_{i+1}^n
    \end{bmatrix}
    \end{split}
    \end{align}

Today we will explore how to formulate a finite differences model that is always stable.

Stability analysis
------------------

Before progressing let us explore when and why the explicit scheme is unstable. Let us assume a temperature profile with a single peak:

.. math::
    :label: eq:imp_dummy2

    T =
    \begin{bmatrix}
    0\\0\\1\\0\\0
    \end{bmatrix}

and let’s further assume that :math:`\beta = 0.2`. Then after one time step the temperature profile looks like 
this:

.. math::
    :label: eq:imp_dummy3

    \begin{align*}
    \begin{split}
    T_2 = \begin{bmatrix}0.2 & 0.6 & 0.2\end{bmatrix}\begin{bmatrix}1\\0\\0\end{bmatrix} = 0.2\\
    T_3 = \begin{bmatrix}0.2 & 0.6 & 0.2\end{bmatrix}\begin{bmatrix}0\\1\\0\end{bmatrix} = 0.6\\
    T_4 = \begin{bmatrix}0.2 & 0.6 & 0.2\end{bmatrix}\begin{bmatrix}0\\0\\1\end{bmatrix} = 0.2
    \end{split}
    \end{align*}


A nice diffusion profile. By further try-and-error (do it for different values for beta!), we find that the scheme is only stable for :math:`\beta<0.5`

Radioactive decay analog
------------------------

A good test problem to explore the stability and accuracy of different finite element schemes is radioactive decay. The change in concentration is controlled by this equation:

.. math::
    :label: eq:radio_1

    \frac{\partial c}{\partial t}= - \lambda c


:eq:`eq:radio_1` states that the rate of decay is controlled by how much of the material is still there. To solve :eq:`eq:radio_1`, we can use our standard finite differences approach:

.. math::
    :label: eq:radio_2

    \frac{c^{n+1} - c^{n}}{\Delta t}= - \lambda c^n


This standard implementation is called the explicit form because you can write the new concentration directly as a function of the old one. The assumption here is that the concentration at the beginning of the time step controls the decay rate.
An alternative formulation, the implicit form, is to assume that the decay rate is controlled by the (unknown) concentration at the end of the time step :math:`c^{n+1}`. The implicit from is therefore:

.. math::
    :label: eq:radio_3

    \frac{c^{n+1} - c^{n}}{\Delta t}= - \lambda c^{n+1}


We will explore the differences between both formulations with a little python script.

FDM notebook
------------

.. toctree::
    :maxdepth: 2

    jupyter/radio_decay.ipynb


Excercise - radioactive decay
-----------------------------
    
Try out the notebook and explore under which conditions the solutions are stable. What's the difference between the explicit and the implicit solution?

Add another method to the script that takes the concentration at the center of the time step:

.. math::
    :label: eq:radio_4

    \frac{c^{n+1} - c^{n}}{\Delta t}= - \lambda c^{n+\frac{1}{2}}

    

Implicit Heat Diffusion
-----------------------

The previous exercise on radioactive decay has shown that the fully implicit method is always stable. We will now rewrite our dike cooling model using the implicit formulation. Here is the implicit finite differences form:

.. math::
    :label: eq:imp_dummy5

    \begin{split}
    \frac{T_i^{n+1} - T_i^n}{\Delta t} = \frac{\kappa}{\Delta x^2} \left( T_{i-1}^{n+1} - 2T_i^{n+1} + T_{i+1}^{n+1} \right)\\
    \begin{bmatrix}
    -\beta & (1+2\beta) & -\beta
    \end{bmatrix}
    \begin{bmatrix}
    T_{i-1}^{n+1}\\
    T_{i}^{n+1}\\
    T_{i+1}^{n+1}
    \end{bmatrix}
    = T_i^{n}
    \end{split}


.. figure:: /_figures/1D_implicit_vs_explicit_stencil.*
   :align: center
   :name: imp_exp_stencil_fig
   :figwidth: 100%

   Different explicit and implicit FD discretization stencils.




It is characteristic for the implicit form that the solution for :math:`T^{n+1}_i` depends on the solution of the neighboring nodes (:math:`i-1` and :math:`i+1`). As a consequence, we cannot simply solve for one node after the other anymore but need to solve a system of equations. In fact, the way forward is to state the problem in matrix form :math:`A\vec{x}=\vec{b}` and solve all equations simultaneously: 

.. math::
    :label: eq:imp_dummy6

    \begin{align}
    \begin{split}
    \begin{bmatrix}
    1 & 0 & 0 & 0 & 0 & 0 & 0\\
    -\beta & (1+2\beta) & -\beta & 0 & 0 & 0 & 0\\
    0 & -\beta & (1+2\beta) & -\beta & 0 & 0 & 0\\
    0 & 0 & -\beta & (1+2\beta) & -\beta & 0 & 0\\
    0 & 0 & 0 & -\beta & (1+2\beta) & -\beta & 0\\
    0 & 0 & 0 & 0 & -\beta & (1+2\beta) & -\beta\\
    0 & 0 & 0 & 0 & 0 & 0 & 1
    \end{bmatrix}
    \begin{bmatrix}
    T_1^{n+1}\\
    T_2^{n+1}\\
    T_3^{n+1}\\
    T_4^{n+1}\\
    T_5^{n+1}\\
    T_6^{n+1}\\
    T_7^{n+1}
    \end{bmatrix}
    =
    \begin{bmatrix}
    T_{top}\\
    T_2^n\\
    T_3^n\\
    T_4^n\\
    T_5^n\\
    T_6^n\\
    T_{bottom}\\
    \end{bmatrix}\\
    \end{split}
    \end{align}


All we need to do is set up the matrix A and then solve for :math:`T^{n+1} = A \backslash T^n`!

FDM notebook
------------

.. toctree::
    :maxdepth: 2

    jupyter/imp_dike_cooling.ipynb


Excercise - implicit dike cooling
---------------------------------
        
    - Get the notebook to work and programm the implicit solution. Explore if you find any stability limits.
    - **Bonus:** Plot the explicit and implicit solutions together. Do you find the same behavior as in the radioactive decay example?