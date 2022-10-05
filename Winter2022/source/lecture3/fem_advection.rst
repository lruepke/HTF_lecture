FEM and advection problems
==========================

In one of the introductory chapters, we had looked into transport problems and how to solve them using FDM. Today we will take this a bit further in order to understand why certain solutions did not work and why finite elements have their own problems for advection. For this purpose, we will again use the standard 1-D steady-state advection diffusion equation:

.. math::
    :label: eq:fem_adv_diff

    c \frac{\partial }{\partial x}T - \frac{\partial}{\partial x}k \frac{\partial T}{\partial x} = 0,

where :math:`c` is velocity, :math:`T` temperature, and :math:`k` thermal conductivity. 

1D Finite Differences Solutions
--------------------------------

Let's have a look at the two different FD solutions that we know alreay: FTCS and upwind. We can write them in condensed form. The FTCS looks like this:

.. math::
    :label: eq:adv_diff_ftcs

    \begin{align}
    \begin{split}
    \frac{c}{\Delta x}
    \begin{bmatrix} 
    -\frac{1}{2} & 0 &  -\frac{1}{2}
    \end{bmatrix} \begin{bmatrix} T_{i-1} \\ T_{i} \\ T_{i+1} \end{bmatrix} -\frac{k}{\Delta x^2} \begin{bmatrix} 
    1 & -2 &  1 \end{bmatrix}\begin{bmatrix} T_{i-1} \\ T_{i} \\ T_{i+1} \end{bmatrix} &= 0 \\
    \begin{bmatrix} 
    -\left(\frac{c}{2\Delta x} + \frac{k}{\Delta x^2} \right) & \frac{2k}{\Delta x^2} &  - \left( - \frac{c}{2 \Delta x} + \frac{k}{\Delta x^2} \right) \end{bmatrix} \begin{bmatrix} T_{i-1} \\ T_{i} \\ T_{i+1} \end{bmatrix}  &= 0 \\
    \end{split}
    \end{align}


Do you notice the asymmetry?

Let's do the same thing for the upwind case:

.. math::
    :label: eq:adv_diff_uw

    \begin{align}
    \begin{split}
    \frac{c}{\Delta x}
    \begin{bmatrix} 
    -1 & 1 &  0
    \end{bmatrix} \begin{bmatrix} T_{i-1} \\ T_{i} \\ T_{i+1} \end{bmatrix} -\frac{k}{\Delta x^2} \begin{bmatrix} 
    1 & -2 &  1 \end{bmatrix}\begin{bmatrix} T_{i-1} \\ T_{i} \\ T_{i+1} \end{bmatrix} &= 0 \\
    \begin{bmatrix} 
    -\left(\frac{c}{\Delta x} + \frac{k}{\Delta x^2} \right) & \left( \frac{c}{\Delta x} + \frac{2k}{\Delta x^2} \right) &  - \frac{k}{\Delta x^2}  \end{bmatrix} \begin{bmatrix} T_{i-1} \\ T_{i} \\ T_{i+1} \end{bmatrix}  &= 0 \\
    \end{split}
    \end{align}


Again asymetric. Let's quickly program those two solutions!


.. toctree::
    :maxdepth: 2

    jupyter/fdm_adv_diff.ipynb


Explore the solution for different pairs of Peclet numbers and c/k ratios. What is the main difference between the upwind and FTCS solution? What kind of numerical diffusion is happening? What's the sign of the numerical diffusion? Is negative diffusion physically possible, or a source of trouble?

1D Finite Element Solutions
--------------------------------
Let’s see how a simple finite element solution looks like. This would be the stiffness matrix for a single element using linear shape functions:

.. math::
    :label: eq:fem_1d_single_element_adv

    \begin{bmatrix}
    \int_{X_e} \left ( N_1 c \frac{\partial N_1}{dx} + \frac{\partial N_1}{\partial x}k\frac{\partial N_1}{\partial x} \right) dx & 
    \int_{X_e} \left (  N_1 c \frac{\partial N_2}{dx} + \frac{\partial N_1}{\partial x}k\frac{\partial N_2}{\partial x} \right) dx \\
    \int_{X_e} \left ( N_2 c \frac{\partial N_1}{dx} + \frac{\partial N_2}{\partial x}k\frac{\partial N_1}{\partial x} \right )dx &
    \int_{X_e} \left ( N_2 c \frac{\partial N_2}{dx} + \frac{\partial N_2}{\partial x}k\frac{\partial N_2}{\partial x} \right )dx
    \end{bmatrix}
    \begin{bmatrix}
    T_1 \\
    T_2
    \end{bmatrix}=0
    
Notice how the sign changed on the diffusion term after partial integration. Again, we get a :math:`2x2` stiffness matrix for every element. The node numbering is again local; node number 1 is the first node of an element k and has the global node number k, while node number 2 is the second node of the element and has the global node number k+1!

We can do the integration using symbolic math. We use again our standard 1D linear interpolation functions; just changed to use :math:`\bar{x}` as our variable describing distance along the element to make the integration easier.

.. math::
    :label: eq:fem_1d_shape_elem_adv
    
    \begin{align}
    \begin{split}
    N_i(x) = 1 - \frac{x-x_i}{x_{i+1} - x_i} = 1 - \frac{x-x_i}{\Delta x} &= 1 - \frac{\bar{x}}{\Delta x}\\
    N_{i+1}(x) = \frac{x-x_i}{x_{i+1} - x_i} = \frac{x-x_i}{\Delta x}& = \frac{\bar{x}}{\Delta x}\\ 
    \end{split}
    \end{align}


Follow the implementation steps in the notebook!

.. toctree::
    :maxdepth: 2

    jupyter/fem_1d_symbolic_adv.ipynb


What we find is that the element stiffness matrix of a single element looks like this:

.. math::
    :label: eq:fem_1d_int_Ae_adv

    Ael = 
    \begin{bmatrix}
    -\frac{c}{2} + \frac{k}{\Delta x} & 
    \frac{c}{2} -\frac{k}{\Delta x} \\
    -\frac{c}{2} -\frac{k}{\Delta x} &
    \frac{c}{2} + \frac{k}{\Delta x}
    \end{bmatrix}

If we assemble it into the global stiffness matrix assuming two elements, we get this equation:

.. math::
    :label: eq:fem_1d_int_Ag_adv

    A = 
    \begin{bmatrix}
    -\frac{c}{2} + \frac{k}{\Delta x}  & \frac{c}{2} -\frac{k}{\Delta x} & 0 \\
    -\frac{c}{2} -\frac{k}{\Delta x} & \frac{2k}{\Delta x} & \frac{c}{2} -\frac{k}{\Delta x} \\
    0 & -\frac{c}{2} -\frac{k}{\Delta x} & \frac{c}{2} + \frac{k}{\Delta x} \\
    \end{bmatrix} 
    \begin{bmatrix}
    T_{1}\\
    T_2\\
    T_{3}\\
    \end{bmatrix}=0

And if we spell this out for the internal node :math:`2`, we get:


.. math::
    :label: eq:fem_1d_int_Ag_adv_node

    \begin{bmatrix}
    -\frac{c}{2} -\frac{k}{\Delta x} & \frac{2k}{\Delta x} & \frac{c}{2} -\frac{k}{\Delta x}
    \end{bmatrix} 
    \begin{bmatrix}
    T_{1}\\
    T_2\\
    T_{3}\\
    \end{bmatrix}=0,

which we can compare to the FTCS solution:

.. math::
    :label: eq:fem_1d_int_Ag_adv_node_ftcs

    \begin{bmatrix} 
    -\left(\frac{c}{2\Delta x} + \frac{k}{\Delta x^2} \right) & \frac{2k}{\Delta x^2} &  - \left( - \frac{c}{2 \Delta x} + \frac{k}{\Delta x^2} \right) \end{bmatrix} \begin{bmatrix} T_{i-1} \\ T_{i} \\ T_{i+1} \end{bmatrix}  &= 0 

 
The two solutions are the same - safe for a constant factor :math:`\Delta x`, which drops out because it also shows up on the RHS (if we had one)! This implies that FE solutions also have an intrinsinc negative diffusion component, which makes them unstable. A popular way to stabilize them is to use the Streamline Upwind Petrov Galerkin (SUPG) method, which modifies the weighting functions so that a little bit of numerical diffusion is added in the direction of flow.