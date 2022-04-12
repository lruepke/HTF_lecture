Introduction: The Finite Differences Method (FDM)
=================================================

Before progressing towards finite element modeling, we will learn about the Finite Difference Method (FDM), which is a somewhat easier method to solve partial differential equations. We will do so by looking at how heat conduction “works”.

Background on conductive heat transport
---------------------------------------

Temperature is one of the key parameters that control and affect the geological processes that we are exploring during this class. Mantle melting, rock rheology, and metamorphism are all functions of temperature. It is therefore essential that we understand how rocks exchange heat and change their temperature.

One fundamental equation in the analysis of heat transfer is Fourier’s law for heat flow:
 
.. math::
    :label: eq:heat_flow
    
    \vec{q} = -k \nabla T = -k
    \begin{bmatrix}
    \frac{\partial T}{\partial x} \\
    \frac{\partial T}{\partial y} \\
    \frac{\partial T}{\partial z}
    \end{bmatrix}

It states that heat flow is directed from high to low temperatures (that’s where the minus sign comes from) and is proportional to the geothermal gradient. The proportionality constant, k, is the thermal conductivity which has units of W/m/K and describes how well a rock transfers heat. k is a typically a complex function of rock type, porosity, and temperature yet is often simplified to a constant value.

In most cases, we are interested in temperature and not heat flow so that we would like to have an equation that describes the temperature evolution in a rock. We can get this by starting from the convervations equation of internal energy.

.. figure:: /_figures/3D_cube_heat_flow.*
   :align: center
   :name: heat_flux_box
   :figwidth: 80%

   Derivation of the energy equation. Change in internal energy is related to changes in heat fluxes into and out of the box (and a lot of other terms (e.g. advection) that are neglected here).

For simple incompressible cases, changes in internal energy can be expressed as changes in temperature times density and specific heat. The change in internal energy with time can now be written as :math:`\rho c_p \frac{\partial}{\partial t} T \Delta x \Delta y \Delta z` (:math:`\rho` is density, :math:`c_p` specific heat, and :math:`T` temperature), has units of (:math:`J/s`), and must be equal to the difference between the heat flow into the box :math:`q_{in} \Delta y \Delta z \left( \frac{W}{m K}\frac{K}{m}m^2 = \frac{J}{s}\right)` and the heat flow out of the box :math:`\left( q_{in} + \frac{\partial q_{in}}{\partial x} \Delta x\right) \Delta y \Delta z` (the y and z directions are done in the same way). With these considerations, we can write a conservation equation for energy in 1D:

.. math::
    :label: eq:1D_heat_flow

    \begin{align}
    \begin{split}
    \rho c_p \frac{\partial T}{\partial t} = -\frac{\partial q_{in}}{\partial x} = \frac{\partial}{\partial x}k\frac{\partial T}{\partial x}\\
    \frac{\partial T}{\partial t} = \frac{k}{\rho c_p} \frac{\partial^2 T}{\partial x^2} = \kappa \frac{\partial^2 T}{\partial x^2}\\
    \frac{\partial T}{\partial t} = \kappa \frac{\partial^2 T}{\partial x^2}   
    \end{split}
    \end{align}

:eq:`eq:1D_heat_flow` is called the heat transfer or heat diffusion equation and is one of the most fundamental equations in Earth Sciences. If the thermal conductivity is constant, we can define a thermal diffusivity, :math:`\kappa`, and write the simpler second form of the equation.

Note that in 3D, the changes in heat flow are expressed as a divergence:

.. math::
    :label: eq:1D_heat_flow_div

    \rho c_p \frac{\partial T}{\partial t} =
    -\nabla \cdot (\vec{q}) =
    \frac{\partial}{\partial x}k\frac{\partial T}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T}{\partial y} + \frac{\partial}{\partial z}k\frac{\partial T}{\partial z}

Or in vector notation:

.. math::
    :label: eq:1D_heat_flow_vec

    \rho c_p \frac{\partial T}{\partial t} =
    -\nabla \cdot \left(\vec{q}\right) =
    \begin{bmatrix}
    \frac{\partial}{\partial x} & \frac{\partial}{\partial y} & \frac{\partial}{\partial z}
    \end{bmatrix}
    k
    \begin{bmatrix}
    \frac{\partial T}{\partial x}\\
    \frac{\partial T}{\partial y}\\
    \frac{\partial T}{\partial z}
    \end{bmatrix}

Finite Differences discretization
---------------------------------

:eq:`eq:1D_heat_flow` is a partial differential equation that describes the evolution of temperature. There are two fundamentally different ways to solve it: 1) analytically or 2) numerically. Analytical solutions have the virtue that they are exact but it is often not possible to find one for complex systems. Numerical solutions are always approximations but can be found also for very complex systems. We will first use one numerical technique called finite differences. To learn how partial differential equations are solved using finite differences, we have to go back to the definition of a derivative:

.. math::
    :label: eq:FD_heat_flow_d1

    \frac{\partial T}{\partial x} =
    \lim_{\Delta x \rightarrow 0} \frac{T(x + \Delta x) - T(x)}{\Delta x}

In our case, the above derivative describes the change in temperature with x (our space coordinate). In numerics, we always deal with discrete functions (as opposed to continuous functions), which only exist at grid points :numref:`dike_setup`. We can therefore translate the above equation into computer readable form:

.. math::
    :label: eq:FD_heat_flow_d2

    \frac{\partial T}{\partial x} \approx
    \frac{T(i + 1) - T(i)}{x(i+1) - x(i)} =
    \frac{T(i + 1) - T(i)}{\Delta x}	

:math:`T(i)` is the temperature at a grid point :math:`i` and :math:`\Delta x` is the grid spacing between two grid points. Using this definition of a derivative, we can now compute the heat flow from the calculated temperature solution:

.. math::
    :label: eq:FD_heat_flow_d3

    q_x = -k \frac{\partial T}{\partial x} \approx -k \left(\frac{T(i+1)-T(i)}{x(i+1)-x(i)}\right)

This form is called the *finite differences* form and is a first step towards solving partial differential equations numerically.

Note that it matters in which direction you count. Usually it makes life much easier if indices and physical coordinates point in the same direction, e.g. x coordinate and index increase to the right.

We have learned how we can compute derivatives numerically. The next step is to solve the heat conduction equation (:eq:`eq:1D_heat_flow`) completely numerically. We are interested in the temperature evolution as a function of time and space :math:`T(x, t)`, which satisfies :eq:`eq:1D_heat_flow`, given an initial temperature distribution. We know already from the heat flow example how to calculate first derivatives (forward differencing):

.. math::
    :label: eq:FD_temperature

    \frac{\partial T}{\partial t} = \frac{T_i^{n+1} - T_i^n}{\Delta t}

The index :math:`n` corresponds to the time step and the index :math:`i` to the grid point (x-coordinate, :numref:`dike_setup`). Next, we need to know how to write second derivatives. A second derivative is just a derivative of a derivate. So we can write (central differencing):

.. math::
    :label: eq:FD_heat_flow_central_difference

    \kappa\frac{\partial^2 T}{\partial x^2} \approx \kappa \frac{\frac{T_{i+1}^n - T_i^n}{\Delta x}-\frac{T_{i}^n - T_{i-1}^n}{\Delta x}}{\Delta x} = \kappa \frac{T_{i+1}^n-2T_{i}^n+T_{i-1}^n}{\Delta x^2}

If we combine equation :eq:`eq:FD_temperature` and :eq:`eq:FD_heat_flow_central_difference` we get:

.. math::
    :label: eq:FD_heat_flow_explicit

    \frac{T_i^{n+1}-T_i^n}{\Delta t} = \kappa \left(\frac{T_{i+1}^n - 2T_i^n + T_{i-1}^n}{{\Delta x^2}}\right)

.. tip::
    Notice how we have *conveniently* used the time index :math:`n` for the temperatures in the spatial derivatives. This results in the explicit form of the final discretized equation. The implicit form, which we will learn about later, would used the unknown new (time index :math:`n+1`) temperatures for the spatial derivatives, which requires solving a system of equations.


The last step is a rearrangement of the discretized equation, so that all known quantities (i.e. temperature at time :math:`n`) are on the right-hand side and the unknown quantities on the left-hand side (properties at :math:`n+1`). This results in:

.. math::
    :label: eq:FD_heat_flow_explicit_solution
    
    T_i^{n+1} = \frac{\kappa \Delta t}{\Delta x^2} \left( T_{i+1}^n - 2 T_i^n + T_{i-1}^n\right) + T_i^n

We have now translated the heat conduction equation :eq:`eq:1D_heat_flow` into a computer readable finite differences form.

Appendix
--------

Taylor-series expansions
^^^^^^^^^^^^^^^^^^^^^^^^

Finite difference approximations can be derived through the use of Taylor-series expansions. Suppose we have a function :math:`f(x)`, which is continuous and differentiable over the range of interest. Let’s also assume we know the value :math:`f(x_0)` and all the derivatives at :math:`x = x_0`. The forward Taylor-series expansion for :math:`f(x_0+\Delta x)` about :math:`x_0` gives

.. math::
    :label: eq:Taylor_series_expansion
    
    f(x_0+\Delta x) =
    f(x_0)+
    \frac{\partial f(x_0)}{\partial x} \Delta x +
    \frac{\partial^2 f(x_0)}{\partial x^2} \frac{(\Delta x)^2}{2!} +
    \frac{\partial^3 f(x_0)}{\partial x^3} \frac{(\Delta x)^3}{3!} +
    \frac{\partial^n f(x_0)}{\partial x^n} \frac{(\Delta x)^n}{n!} +
    O(\Delta x)^{n+1}

We can compute the first derivative by rearranging equation \ref{eqs:Taylor_series_expansion}

.. math::
    :label: eq:Taylor_series_expansion_rearranged

    \frac{\partial f(x_0)}{\partial x} =
    \frac{f(x_0+\Delta x)− f(x_0)}{\Delta x} −
    \frac{\partial^2 f(x_0)}{\partial x^2} \frac{(\Delta x)}{2!} −
    \frac{\partial^3 f(x_0)}{\partial x^3} \frac{(\Delta x)^2}{3!} ...

This can also be written in discretized notation as:

.. math::
    :label: eq:Taylor_series_expansion_rearranged_2

    \frac{\partial f(x_i)}{\partial x} = \frac{f_{i+1}−f_i}{\Delta x} + O(\Delta x)

here :math:`O(\Delta x)` is called the *truncation error*, which means that if the distance :math:`\Delta x` is made smaller and smaller, the (numerical approximation) error decreases with :math:`\Delta x`. This derivative is also called first order accurate.

We can also expand the Taylor series backward

.. math::
    :label: eq:Tayler_series_expansion_backward

    f(x_0−\Delta x) =
    f(x_0)− \frac{\partial f(x_0)}{\partial x} \Delta x+
    \frac{\partial^2 f(x_0)}{\partial x^2} \frac{(\Delta x)^2}{2!} −
    \frac{\partial^3 f(x_0)}{\partial x^3} \frac{(\Delta x)^3}{3!} + \cdots

In this case, the first (backward) derivative can be written as

.. math::
    :label: eq:Tayler_series_expansion_backward_fist_derivativ

    \frac{\partial f (x0)}{\partial x} =
    \frac{f(x_0)− f(x_0−\Delta x)}{\Delta x} +
    \frac{\partial^2 f(x_0)}{\partial x^2} \frac{(\Delta x)}{2!} −
    \frac{\partial^3 f(x_0)}{\partial x^3} \frac{(\Delta x)^2}{3!} \cdots\\
    \frac{\partial f(x_i)}{\partial x} = \frac{f_i− f_{i−1}}{\Delta x} + O(\Delta x)

By adding :eq:`eq:Taylor_series_expansion_rearranged` and :eq:`eq:Tayler_series_expansion_backward_fist_derivativ` and dividing by two, a second order accurate first order derivative is obtained

.. math::
    :label: eq:Tayler_series_expansion_dummy1

    \frac{\partial f(x_i)}{\partial x} = \frac{f_{i+1} − f_{i−1}}{2\Delta x} + O(\Delta x)^2

By adding equations :eq:`eq:Taylor_series_expansion_rearranged` and :eq:`eq:Tayler_series_expansion_backward` an approximation of the second derivative is obtained

.. math::
    :label: eq:Tayler_series_expansion_dummy2

    \frac{\partial f^2(x_i)}{\partial x^2} = \frac{f_{i+1}−2 f_i+ f_{i−1}}{(\Delta x)^2} +O(\Delta x)^2

With this approach we can basically derive all possible finite difference approximations. A different way to derive the second derivative is by computing the first derivative at :math:`i+\frac{1}{2}` and at :math:`i-\frac{1}{2}` and computing the second derivative at :math:`i` by using those two first derivatives:

.. math::
    :label: eq:Tayler_series_expansion_dummy3
    
    \begin{align}
    \begin{split}
    \frac{\partial f(x_{i+1/2})}{\partial x} = \frac{f_{i+1}−f_i}{x_{i+1}−x_i}\\
    \frac{\partial f(x_{i−1/2})}{\partial x} = \frac{f_i−f_{i−1}}{x_i−x_{i−1}}\\
    \frac{\partial f^2(x_i)}{\partial x^2} = \frac{\frac{\partial f(x_{i+1/2})}{\partial x} − \frac{\partial f(x_{i−1/2})}{\partial x}}{x_{i+1/2}−x_{i−1/2}} =
    \frac{ \frac{f_{i+1}−f_i}{x_{i+1}−x_i} − \frac{f_i−f_{i−1}}{x_i−x_{i−1}}}{0.5(x_{i+1}−x_{i−1})}
    \end{split}
    \end{align}

Similarly we can derive higher order derivatives. Note that the highest order derivative that usually occurs in geodynamics is the :math:`4^th`-order derivative.


Finite difference approximations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following equations are common finite difference approximations of derivatives. If you, in the future, need to write a finite difference approximation, come back here.

Left-sided first derivative, first order

.. math::
    :label: eq:Tayler_series_expansion_dummy4
    
    \left| \frac{\partial u}{\partial x}\right|_{i−1/2} = \frac{u_i−u_{i−1}}{\Delta x} +O(\Delta x)

Right-sided first derivative, first order

.. math::
    :label: eq:Tayler_series_expansion_dummy5
    
    \left| \frac{\partial u}{\partial x}\right|_{i+1/2} = \frac{u_{i+1}−u_i}{\Delta x} +O(\Delta x)

Central first derivative, second order

.. math::
    :label: eq:Tayler_series_expansion_dummy6
    
    \left| \frac{\partial u}{\partial x}\right|_i = \frac{u_{i+1}−u_{i−1}}{2\Delta x} +O(\Delta x)^2

Central first derivative, fourth order

.. math::
    :label: eq:Tayler_series_expansion_dummy7
    
    \left| \frac{\partial u}{\partial x}\right|_i = \frac{−u_{i+2}+8u_{i+1}−8u_{i−1}+u_{i−2}}{12\Delta x} +O(\Delta x)^4

Central second derivative, second order

.. math::
    :label: eq:Tayler_series_expansion_dummy8
    
    \left| \frac{\partial^2 u}{\partial x^2}\right|_i = \frac{u_{i+1}−2u_i+u_{i−1}}{\Delta x^2} +O(\Delta x)^2

Central second derivative, fourth order

.. math::
    :label: eq:Tayler_series_expansion_dummy9
    
    \left|\frac{\partial^2 u}{\partial x^2} \right|_i = \frac{−u_{i+2}+16u_{i+1}−30u_i+16u_{i−1}−u_{i−2}}{12\Delta x^2} +O(\Delta x)^4

Central third derivative, second order

.. math::
    :label: eq:Tayler_series_expansion_dummy10
    
    \left|\frac{\partial^3 u}{\partial x^3}\right|_i = \frac{u_{i+2}−2u_{i+1}+2u_{i−1}−u_{i−2}}{2\Delta x^3} +O(\Delta x)^2

Central third derivative, fourth order

.. math::
    :label: eq:Tayler_series_expansion_dummy11
    
    \left| \frac{\partial^3 u}{\partial x^3}\right|_i = \frac{−u_{i+3}+8u_{i+2}−13u_{i+1}+13u_{i−1}−8u_{i−2}+u_{i−3}}{8\Delta x^3} +O(\Delta x)^4

Central fourth derivative

.. math::
    :label: eq:Tayler_series_expansion_dummy12
    
    \left| \frac{\partial^4 u}{\partial x^4}\right|_i =\frac{u_{i+2}−4u_{i+1}+6u_i−4u_{i−1}+u_{i−2}}{\Delta x^4} +O(\Delta x)2

Note that the higher the order of the finite difference scheme, the more adjacent points are required. It is also important to note that derivatives of the following form

.. math::
    :label: eq:Tayler_series_expansion_dummy13
    
    \frac{\partial}{\partial x}\left(k\frac{\partial u}{\partial x}\right)

should be formed as follows

.. math::
    :label: eq:Tayler_series_expansion_dummy14
    
    \left| \frac{\partial}{\partial x}\left(k\frac{\partial u}{\partial x} \right)\right|_i =
    \frac{k_{i+1/2} \frac{u_{i+1}−u_i}{\Delta x} − k_{i−1/2} \frac{u_i−u_{i−1}}{\Delta x}}{\Delta x} +O(\Delta x)^2


