Transport problems with FDM
=================================================

Before progressing towards finite elements, we will look into transport problems using FDM. Advection problems are described by so-called hyperbolic partial differential equations. The general form of advection/transport equations looks like this:

.. math::
    :label: eq:advection_eqn
    
    \frac{\partial C}{\partial t} = -\vec{u}\cdot\nabla C

or written in one-dimension:

.. math::
    :label: eq:advection_eqn_1d
    
    \frac{\partial C}{\partial t} = -u_x \frac{\partial C}{\partial x}

These equations describe the non-diffusive transport of a concentration :math:`C` with the flow field :math:`\vec{u}`, like, for example, in solute transport. If the solute also diffused, there would also be a diffusion term on the right-hand side. 

While these equations look quite "friendly" and a bit simpler than the diffusion equations we have looked at so far, they can actually be a challenge to solve (well). Let's explore a two possible FDM schemes!


FDM advection schemes
---------------------------------------

Forward in time, central in space (FTCS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
One possible FDM discretization is to use the so-called *forward in time, central in space (FTCS)* technique, which uses forward differencing for the time derivative (what we had also done in the diffusion case) and central differencing for the gradient term.


.. math::
    :label: eq:advection_ftcs
    
    \frac{C_i^{n+1}-C_i^n}{\Delta t}= -u_{x,i} \frac{C_{i+1}^{n}-C_{i-1}^n}{2\Delta x}


Upwind
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Another option is to use what is called an upwind scheme, in which the gradient on the right-hand side is evaluated according to the direction of flow:

.. math::
    :label: eq:advection_upwind
    
    \frac{C_i^{n+1}-C_i^n}{\Delta t}= -u_{x,i}  \begin{pmatrix} \frac{C_{i}^{n}-C_{i-1}^n}{\Delta x} \text{ , if $u_x \geq 0$} \\ \frac{C_{i+1}^{n}-C_{i}^n}{\Delta x} \text{ , if $u_x < 0$} \end{pmatrix} 


There are many more schemes including much better ones (see e.g. Marc Spiegelmann's script on numerical methods  `<https://earth.usc.edu/~becker/teaching/557/reading/spiegelman_mmm.pdf>`_). We here look into these two basic types because we will come back to them when looking into finite elements for advection problems.

Excercise
^^^^^^^^^^

Implement the different schemes using this jupyter notebook as a starting point. Notice how the FTCS scheme is always unstable and how the upwind scheme is very diffusive.

.. toctree::
    :maxdepth: 2

    jupyter/fdm_advection_schemes.ipynb