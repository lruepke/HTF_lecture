Introduction Finite Element Method (FEM)
=========================================

We will again use the 1-D steady state heat diffusion equation as an example and solve it using FEM. Formally speaking the steady-state diffusion equation is an **elliptic** PDE. Elliptic equations describe a large variety of steady state problems in earth sciences including diffusion, plate flexure, and incompressible flows.

.. admonition:: Types of PDEs

    Besides **elliptic** PDEs, there are, for example, also **parabolic** PDEs, which include the transient term, and **hyperbolic** PDEs, which describe e.g. wave propagation. We will learn about those a bit later in the course.

Governing equation - strong form
--------------------------------

Let's go on and solve the elliptic 1-D heat diffusion equation with the Finite Element Method (FEM):

.. math::
    :label: eq:fem_1d_strong

    \frac{\partial}{\partial x}k\frac{\partial T_{ex}}{\partial x}=0. 

This is called the strong form of the governing equation and the function :math:`T_{ex}(x)`is the exact solution to it. We now follow our MWR recipe and replace the continuous function :math:`T_{ex}(x)` by a (simpler) approximation:


Approximate solution
--------------------
.. math::
    :label: eq:fem_aprox_funcion

    T_{ex} \cong \tilde{T}= \sum_{j=1}^{n} N_j T_j = N_j T_j. 

In FEM we divide our computational domain into simple geometeric entities, the finite elements, and those finite elements have grid points. For example, in 2-D we can divide the domain in triangles with three nodes and in 1-D in simple lines with two nodes each. Think of our simple FD descritization, if we were to use the same descritization in FEM, the intervals between grid points would be our finite elements.

.. figure:: /_figures/shapeFunction_Linear_1D.svg
    :name: fig:shapeFunc:1D:linear
    :align: center

    Example of 1D linear finite element shape functions

These linear independent functions :math:`N`, called shape functions in FEM, are associated with the chosen discretization (:numref:`fig:shapeFunc:1D:linear`). Each node is associated with one shape function and that shape function has the value :math:`1` at that node and is zero at every other node. In our simple 1-D case, the shape functions vary linearly from one to zero over the adjacent elements of each node. With these properties, the sum over all all shape functions is always one, and the coefficients of those shape functions happen to be the unknown variable values at the nodes, the nodal temperatures in our case.

Let's plug the approximated solution into the strong form of our equation:

.. math::
    :label: eq:fem_1D_R

    R = \frac{\partial}{\partial x}k\frac{\partial N_j T_j }{\partial x}\neq 0. 

Note how we use the Einstein convention and sum over repeated indices. We can now write the MWR form of the equation:

.. math::
    :label: eq:fem_mwr
    
    \int_XW_i\frac{\partial}{\partial x}k\frac{\partial N_j T_j }{\partial x}dx=0\ \ \ \ \ \ \ i=1,2,...,n


Galerkin method
----------------
In finite elements, the most frequently used weighting method is the Galerkin method that we already know from the previous session. We therefore use the interpolation functions also as weighting functions:

.. math::
    :label: eq:fem_mwr_2
    
    \int_XN_i\frac{\partial}{\partial x}k\frac{\partial N_j T_j }{\partial x}dx=0\ \ \ \ \ \ \ i=1,2,...,n

The next step is to reduce the order of differentiation by using partial integration. Remember, integration by parts looks like this:

.. math::
    :label: eq:partial_int

    \int_{\Omega} fg' d\Omega = -\int_{\Omega} f'g d\Omega + \oint_{\Gamma} fg d\Gamma


Weak form
----------
Applying this to :eq:`eq:partial_int` results in the **weak form** of our governing equation:

.. math::
    :label: eq:fem_1d_weak

    \begin{align}
    \begin{split}
    \int_XN_i\frac{\partial}{\partial x}k\frac{\partial N_j T_j }{\partial x}dx=0\ \ \ \ \ \ \ i=1,2,...,n\\
    \Rightarrow \\
    -\int_X \frac{\partial N_i}{\partial x}k\frac{\partial N_j T_j }{\partial x}dx + \oint_{\Gamma} N_ik\frac{\partial N_j T_j }{\partial x}d\Gamma=0\ \ \ \ \ \ \ i=1,2,...,n\\
    \Rightarrow \\
    \int_X \frac{\partial N_i}{\partial x}k\frac{\partial N_j T_j }{\partial x}dx + \oint_{\Gamma}\vec{q}\vec{n}d\Gamma=0\ \ \ \ \ \ \ i=1,2,...,n\\
    \end{split}
    \end{align}

Close inspection of the boundary integral reveals that it is the heat flow through the boundaries of the modeling domain. In 1-D these are simply numbers of heat flow at the sides and in 2-D/3-D heat fluxes through the sides of the modeling domain. For the moment we will neglect it and only continue with the solution inside our domain:

.. math::
    :label: eq:fem_1d_weak_simple

    \int_X \frac{\partial N_i}{\partial x}k\frac{\partial N_j T_j }{\partial x}dx =0\ \ \ \ \ \ \ i=1,2,...,n

How can we solve this using the FEM? The basic idea is to split the integrals into subdomains and the solution will be the sum of the subdomains. In FEM, those subdomains are the finite elements:
    
.. math::
    :label: eq:fem_1d_weak_simple_2

    \begin{align}
    \begin{split}
    \int_X \frac{\partial N_i}{\partial x}k\frac{\partial N_j T_j }{\partial x}dx =0\ \ \ \ \ \ \ i=1,2,...,n\\
    \Rightarrow \\
    \int_X \frac{\partial N_i}{\partial x}k\frac{\partial N_j T_j }{\partial x}dx =  \sum_{Elements} \int_{X_e} \frac{\partial N_i}{\partial x}k\frac{\partial N_j T_j }{\partial x}dx = 0\ \ \ \ \ \ \ i=1,2,...,n\\
    \end{split}
    \end{align}


While mathematically we are always allowed to split an integral into a sum of integrals, in FEM this is particularly useful. The catch is that the shape function associated with a node is only non-zero in the elements connected to that node. This in turn means that in :eq:`eq:fem_1d_weak_simple_2`, we can do the integration of each element completely independent from all other elements!

In the following sessions, we will go through all the steps involved in solving :eq:`eq:fem_1d_weak_simple_2`. We will first do this in a simple 1-D case before moving on to the general 2-D.

FEM implementation
==================

Element stiffness matrix :math:`A_{el}`
---------------------------------------
Let’s look at a single 1-D element. If we use linear shape functions, one 1-D element has two nodes – and only the shape functions of those two nodes will be non-zero. The so-called connectivity, which connects elements and nodes (which nodes belong to which element) is easy in 1-D: element number k will have the nodes i=k and i=k+1 (Element 1 has nodes 1 and 2). This means that for each element we will get two equations:

.. math::
    :label: eq:fem_1d_single_element

    \begin{bmatrix}
    \int_{X_e} \frac{\partial N_1}{\partial x}k\frac{\partial N_1}{\partial x}dx & 
    \int_{X_e} \frac{\partial N_1}{\partial x}k\frac{\partial N_2}{\partial x}dx \\
    \int_{X_e} \frac{\partial N_2}{\partial x}k\frac{\partial N_1}{\partial x}dx &
    \int_{X_e} \frac{\partial N_2}{\partial x}k\frac{\partial N_2}{\partial x}dx
    \end{bmatrix}
    \begin{bmatrix}
    T_1 \\
    T_2
    \end{bmatrix}=0
    
We get a :math:`2x2` so-called stiffness matrix for every element. In the end we will have to sum the contributions from every element into a global stiffness matrix, which will again be :math:`[n x n]`.  Note: The node numbering in :eq:`eq:fem_1d_single_element` is local! I.e. node number 1 is the first node of an element k and has the global node number k, while node number 2 is the second node of the element and has the global node number k+1!

We have two linear shape functions for the points i and i+1 of each element (see :numref:`fig:shapeFunc:1D:linear`):

.. math::
    :label: eq:fem_1d_shape_elem
    
    \begin{align}
    \begin{split}
    N_i(x) = 1 - \frac{x-x_i}{x_{i+1} - x_i} &= 1 - \frac{x-x_i}{\Delta x}\\
    N_{i+1}(x) = \frac{x-x_i}{x_{i+1} - x_i} &= \frac{x-x_i}{\Delta x}\\ 
    \end{split}
    \end{align}


Excercise: derive :math:`A_{el}`
----------------------------------
.. toctree::
    :maxdepth: 2

    jupyter/fem_1d_symbolic.ipynb


Matrix assembly
---------------
The above exercise results in the integrated element stiffness matrix of steady-state diffusion:

.. math::
    :label: eq:fem_1d_int_Ae

    Ael = 
    \begin{bmatrix}
    \frac{k}{\Delta x} & 
    -\frac{k}{\Delta x} \\
    -\frac{k}{\Delta x} &
    \frac{k}{\Delta x}
    \end{bmatrix}


The global stiffness matrix :math:`A` is the sum of all element stiffness matrices; we simply have to add them all together using the connectivity information, i.e. which nodes belong to which elements (something that we call EL2NOD in the codes).

If we assume three linear elements, this looks like this

.. math::
    :label: eq:fem_1d_matrix_assembly

    A = 
    \begin{bmatrix}
    \frac{k}{\Delta x} & -\frac{k}{\Delta x} & 0 & 0\\
    -\frac{k}{\Delta x} & \frac{k}{\Delta x} & 0 & 0\\
    0 & 0 & 0 & 0\\
    0 & 0 & 0 & 0\\
    \end{bmatrix} + 
    \begin{bmatrix}
    0 & 0 & 0 & 0\\
    0 & \frac{k}{\Delta x} & -\frac{k}{\Delta x} & 0 \\
    0 & -\frac{k}{\Delta x} & \frac{k}{\Delta x} & 0\\
    0 & 0 & 0 & 0\\
    \end{bmatrix}+ 
    \begin{bmatrix}
    0 & 0 & 0 & 0\\
    0 & 0 & 0 & 0\\
    0 & 0 & \frac{k}{\Delta x} & -\frac{k}{\Delta x} \\
    0 & 0 & -\frac{k}{\Delta x} & \frac{k}{\Delta x} \\
    \end{bmatrix}


Implementation
---------------

Let's put everything together and write our first finite element code! We will do this by programming a FEM solution to the steady-state heat diffusion equation but this time with a source term.

.. math::
    :label: eq:fem_1d_strong_ex

    \frac{\partial}{\partial x}k\frac{\partial T_{ex}}{\partial x}+Q=0.
    
The weak form will then look like this (signs get flipped during partial integration):

.. math::
    :label: eq:fem_1d_weak_simple_ex1

    \begin{align}
    \begin{split}
    \int_X \frac{\partial N_i}{\partial x}k\frac{\partial N_j  }{\partial x}T_j - N_iQdx =0\ \ \ \ \ \ \ i=1,2,...,n\\
    \Rightarrow \\
    \int_X \frac{\partial N_i}{\partial x}k\frac{\partial N_j  }{\partial x} T_j - N_i Qd x =  \sum_{Elements} \int_{X_e} \frac{\partial N_i}{\partial x}k\frac{\partial N_j  }{\partial x} T_j - N_i Q dx = 0\ \ \ \ \ \ \ i=1,2,...,n\\
    \end{split}
    \end{align}

After another symbolic integration (try it) of the source term,we get the following element stiffness matrix and element system of equations:

.. math::
    :label: eq:fem_1d_weak_simple_ex2

    Ael = 
    \begin{bmatrix}
    \frac{k}{\Delta x} & 
    -\frac{k}{\Delta x} \\
    -\frac{k}{\Delta x} &
    \frac{k}{\Delta x}
    \end{bmatrix} 
    \begin{bmatrix}
    T_i \\
    T_{i+1}
    \end{bmatrix} = 
    \begin{bmatrix}
    \frac{\Delta x}{2}Q_{el} \\
    \frac{\Delta x}{2}Q_{el}
    \end{bmatrix}

Where the first matrix is the element stiffness matrix, the second vector are the unknown temperatures of the element, and the right-hand side has the integration source term :math:`Q_{el}` which here is constant value (per element).


Let's put this into python!

.. toctree::
    :maxdepth: 2

    jupyter/fem_1d_numerical.ipynb
