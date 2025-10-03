.. include:: /include.rst_

Porous Flow - Submarine Hydrothermal Systems
============================================

What is the permeability of seafloor?
--------------------------------------
We have just learned how to compute the effective permeability of a rock sample from direct numerical simulations of flow on the pore scale. It is one way of getting the permeability of a rock sample. However, often the effective permeability of a rock sample is not sufficient to describe the flow behavior in a larger geological context, such the hydrothermal circulation pattern that feed black smoker vent sites at the seafloor.

But what is the permeability of the seafloor? The answer is, of course, that depends... It depends on the rock type, the scales we are interested in, and and many other factors. Popular methods of estimating the permeability of the seafloor include:

- laboratory measurements on rock samples,
- field measurements of in-situ permeability using packer tests on IODP :cite:`Fisher1998_permeability`
- seismic methods that infer porosity and relate it to permeability using empirical relations :cite:`Marjanovic2019_permeability`.
- and even methods that use the ocean tides to estimate the bulk permeability of the oceanic crust from phase lags between the tidal forcing and the response of the hydrothermal system :cite:`Barreyre2018_permeability`

Do you remember what we said about scales and homogenization in porous media? The permeability of a rock sample is not a constant, but depends on the scale we are looking at and all these methods look at different scales. A rock sample may have a permeability of 10^-12 m^2, but the bulk permeability of the seafloor may be 10^-14 m^2 or even lower.

This is why there are endless debates about what the "correct" permeability of the seafloor is. The truth is, there is no single correct answer. It depends on the context and the scales we are interested in.


Theory: Hydrothermal convection
-------------------------------------

Let's proceed and assume that we have agreed on a permeability structure and now want to explore how hydrothermal convection actually looks like. During this lecture we will study single-phase hydrothermal flow in submarine hydrothermal systems. The respective solver is named :code:`HTFoam`. The hydrothermal fluid flow is governed by Darcy's law (Eqn. :eq:`eq:darcy_l04`), mass continuity (Eqn. :eq:`eq:conti_l04`) and energy conservation (Eqn. :eq:`eq:temperature`) equations shown below, 

.. math::
    :label: eq:darcy_l04 

    \vec{U} = - \frac{k}{\mu_f} (\nabla p -\rho \vec{g})
      
.. math::
    :label: eq:conti_l04
    
    \varepsilon \frac{\partial \rho_f}{\partial t} + \nabla \cdot (\vec{U} \rho_f) = 0

.. math::
    :label: eq:rho_to_pt

    \frac{\partial \rho}{\partial t} =   \frac{\partial \rho}{\partial P}\frac{\partial P}{\partial t}\rvert_T +  \frac{\partial \rho}{\partial T}\frac{\partial T}{\partial t}\rvert_P

.. math::
    :label: eq:beta

    \beta = \frac{1}{\rho}\frac{\partial \rho}{\partial P}

.. math::
    :label: eq:alpha
    
    \alpha = -\frac{1}{\rho}\frac{\partial \rho}{\partial T}

.. math::
    :label: eq:pressure
    
    \varepsilon \rho_f \left( \beta_f \frac{\partial p}{\partial t} - \alpha_f \frac{\partial T}{\partial t} \right) = \nabla \cdot \left( \rho_f \frac{k}{\mu_f} (\nabla p - \rho_f \vec{g}) \right)

.. math::
    :label: eq:temperature
    
    (\varepsilon \rho_f C_{pf} + (1-\varepsilon)\rho_r C_{pr})\frac{\partial T}{\partial t} = \nabla \cdot (\lambda_r \nabla T) - \rho_f C_{pf} \vec{U}\cdot \nabla T + \frac{\mu_f}{k} \parallel \vec{U} \parallel ^2 - \left( \frac{\partial ln \rho_f}{\partial ln T} \right)_p \frac{Dp}{Dt}

where the pressure equation :eq:`eq:pressure` is derived from continuity equation :eq:`eq:conti_l04` and Darcy's law :eq:`eq:darcy_l04`.

Implementation
--------------
The details of the OpenFoam implementation can be found in the |foam| documentation. Here we only show a brief summary. :numref:`fig:htf_solution` shows how the energy equation is solved within the OpenFoam framework.


.. figure:: /_figures/solution_algorithm.* 
   :align: center
   :name: fig:htf_solution

   Implementation of the energy conservation equation.

Equation-of-state
-----------------
The fluid properties like density, viscosity, specific heat are determined from the equation-of-state of pure water. :numref:`fig:phase_diagram` shows the phase diagram of pure water. At sub-critical conditions (P< 22 MPa), the boiling curve divides the regions of liquid water and water vapor. At super-critical conditions, there is a gradual transition from a liquid-like to a vapor-like fluid phase. |foam| is a single phaes code and can only be used in regions, where a single fluid phase is present, i.e. under pure liquid water, water vapor, or supercritical conditions; boiling cannot be resolved. As we will find out later, the thermodynamic properties of water have first order control on flow dynamics and upflow temperatures in submarine hydrothermal systems. 


.. figure:: /_figures/PhaseDiagram.*
   :align: center
   :name: fig:phase_diagram

   Phase diagram of pure water.
