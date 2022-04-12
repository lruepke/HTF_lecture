Example 3: Periodic variations in seafloor temperature
======================================================

Flux boundary conditions
------------------------

So far, we have always assumed that the boundary condition for the temperature equation was a fixed temperature. This condition is also called a *Dirichlet boundary condition*. We can, however, also assume a case where the boundary has a constant gradient (e.g. heat flux). This is called a *Neumann boundary conditions* and an example would be a given conductive heat flux from below. A flux boundary condition at the bottom could be written like this:

.. math::
    :label: eq:sf_dummy1

    \frac{\partial T(y=-L,t)}{\partial y}=c_1


where :math:`c_1` is the specified gradient. We could also write it as a heat flow boundary conditions

.. math::
    :label: eq:sf_dummy2

    \frac{−q_{base}}{k}= \frac{\partial T(y = −L,t)}{\partial y}

We can program these conditions by using a forward or a backward finite difference expression. However, this is not so good, since these finite difference approximations are only first order accurate in space. Moreover they would yield the first derivative at the location :math:`T_{\frac{1}{2}}` or at :math:`T_{ny-\frac{1}{2}}`, and not at :math:`T_1` and :math:`T_y`. A better way to incorporate a flux boundary conditions is therefore to use a central finite difference approximation, which is given (at :math:`i = 1`) by:

.. math::
    :label: eq:flux_BC

    \frac{T_2^{n+1} − T_0^{n+1}}{2\Delta y} = −\frac{q_{base}}{k}

The problem is that the expression above involves a point that is not part of the numerical grid (:math:`T_0`). A way around this can be found by noting that the equation for the center nodes is given by:

.. math::
    :label: eq:centered_node
    
    \frac{T_i^{n+1} − T_i^n}{\Delta t} = \kappa \frac{T_{i+1}^{n+1} − 2T_i^{n+1} + T_{i−1}^{n+1}}{\Delta y^2}

Writing this expression for the first node gives

.. math::
    :label: eq:centered_node_first

    \frac{T_1^{n+1} − T_1^n}{\Delta t} = \kappa \frac{T_2^{n+1} − 2T_1^{n+1} + T_0^{n+1}}{\Delta y^2}


An explicit expression for :math:`T_0` can be obtained from :eq:`eq:flux_BC`:

.. math::
    :label: eq:centered_node_first_explicit
    
    T_0^{n+1} = T_2^{n+1} + \frac{2 \Delta y q_{base}}{k}


Substituting :eq:`eq:centered_node_first_explicit` into :eq:`eq:centered_node_first` gives:

.. math::
    :label: eq:sf_dummy3
    
    \frac{T_1^{n+1} − T_1^n}{\Delta t} = \kappa \frac{2 T_2^{n+1} − 2 T_1^{n+1} + 2\Delta y \frac{q_{base}}{k}}{\Delta y^2}


Again we can rearrange this equation to bring known terms to the right-hand side:

.. math::
    :label: eq:flux_BC_implicit_formulation

    \begin{align}
    \begin{split}
    (1+2\beta) T_1^{n+1} − 2\beta T_2^{n+1} = T_1^n + 2\beta \Delta y \frac{q_{base}}{k}\\
    \beta = \frac{\kappa \Delta t}{\Delta y^2}
    \end{split}
    \end{align}


This equation only involves grid points that are part of the computational grid and equation :eq:`eq:flux_BC_implicit_formulation` can be incorporated into the matrix :math:`A` and the right-hand side :math:`b`.

Example: seafloor temperature variations
----------------------------------------

Oceanic heat flow measurements provide important insights into cooling and alteration processes of oceanic plates. Such measurements are usually done with devices that measure heat flow within the first few meters of sediment. A natural question to ask is to which degree such measurements may be perturbed by seasonal variations in bottom water temperatures. Let’s set up a simple model for this.

Assume a vertical modeling domain of :math:`30m`. At the top of the domain a sinusoidal change in surface temperature (:math:`\pm 2°C` around :math:`4°C`) over a one year period is applied and at the bottom a constant heat flow of :math:`60 \frac{mW}{m^2}` is assumed. The sediments have a constant diffusivity of :math:`1e^{-6} m^2/s` and a thermal conductivity of :math:`1.5 W/m/K`.

FDM notebook
------------

.. toctree::
    :maxdepth: 2

    jupyter/seafloor_temp_variations.ipynb


Excercise - periodic changes in seafloor temperature
----------------------------------------------------
        
Explore the solution and understand what controls the depth to which periodic varations in seafloor temperature are "felt" inside the sediments

    - To which depth to seasonal, annual, and decadal variations propagate into the sediments?
    - Is there a phase shift between surface heat flow and temperature variations?