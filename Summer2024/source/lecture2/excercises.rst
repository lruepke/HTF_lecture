Excercises
==========

Let us program some of the above MWR formulations. Here is a jupyter notebook that we can use as a starting point.

Notebook - excercises
------------------------------------------
        
.. toctree::
    :maxdepth: 2
    
    jupyter/MWR_excercise.ipynb

Excercises:

    - Get the different methods to work and explore how good the MWR solutions are.
    - What happens if you increase the velocity and make the problem more advection dominated?
    - What do you think is the problem with advectino dominated problems?

Higher-order approximations
----------------------------

What we have seen is that for diffusion dominated problems, the above approximation gives reasonable results. For advection dominated problems, a simple polynomial function can no longer approximate the solution in a satisfactory way. One possible solution is to use a higher polynomial as the approximating function. Let’s try this:

.. math::
    :label: eq:mwr_example_ho_1
    
    \tilde{u}(x)=a_0 + a_1 x + a_2 x^2 + a_3 x^3


The boundary conditions can now be written as:

.. math::
    :label: eq:mwr_example_ho_2

    \begin{align}
    \begin{split}
    \tilde{u}(0) &= 0=a_0\\
    \tilde{u}(1) &= 1=0+a_1 + a_2 + a_3\\
    0 &= a_1 + a_2 + a_3 -1
    \end{split}
    \end{align}  

The first and second derivatives become:

.. math::
    :label: eq:mwr_example_ho_3

    \begin{align}
    \begin{split}
    \frac{\partial \tilde{u}}{\partial x} &= a_1+2a_2x + 3a_3x^2\\
    \frac{\partial^2 \tilde{u}}{\partial x^2} &= 2a_2 + 6a_3x
    \end{split}
    \end{align}   

So the residual is

    .. math::
        :label: eq:mwr_example_ho_4
    
        R(x)=c\left ( a_1 + 2a_2x+3a_3x^2 \right ) - K(2a_2+6a_3x)

We will only use the collocation method with two collocation points to determine the two unknown coefficients. The result is a system of equations for the four coefficients :math:`ai`:

.. math::
    :label: eq:mwr_example_ho_5

    \begin{align}
    \begin{split}
    a_0 &=0\\
    c \left( a_1+2a_2 \left(\frac{1}{3}\right) +3a_3\left(\frac{1}{3}\right)^2 \right)-K\left(2a_2+6a_3\left(\frac{1}{3}\right)\right)&=0=R\left(\frac{1}{3}\right)\\
    c \left( a_1+2a_2 \left(\frac{2}{3}\right) +3a_3\left(\frac{2}{3}\right)^2 \right)-K\left(2a_2+6a_3\left(\frac{2}{3}\right)\right)&=0=R\left(\frac{2}{3}\right)\\
    a_1 + a_2 + a_3 -1 &=0
    \end{split}
    \end{align}   

We can solve this system of equations with the Symbolic Math Toolbox of Matlab or again via sympy using python. Try to do that!

Notebook - derivation of higher-order constants
-----------------------------------------------
        
.. toctree::
    :maxdepth: 2

    jupyter/MWR_higher_order.ipynb

Have a close look at the above notebook and make sure you understand how symbolic systems of equations are solved using sympy! Also check if we are now able to approximate the solution for advection-dominated problems.

Outlook
--------

The above considerations have shown how we can find approximate solutions to partial differential equations using MWR. If you remember the characterization above: *“...finite element methods partition the problem domain into subregions and produce a simple function in each subregion that approximates the solution”*, we can already see that finite elements do something very similar. In essence we will partition the domain into many elements, approximate the solution over every element using the Galerkin method and very simple basis functions, and will finally obtain the global solution by summing up the contributions of each element. That is what we will do over the next weeks!