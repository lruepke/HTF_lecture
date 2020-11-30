.. include:: /include.rst_

Porous Flow - Submarine Hydrothermal Systems
============================================

Theoretical background
-------------------------------------

In the last lecture, we had made our first Navier-Stokes simulations, which resolved the full dynamics of fluid flow. When studying porous media, flow is often approximated by Darcy's law that states that flow is proportional to the pressure gradient. The invovled constants are viscosity and permeability. Extensive theoretical work helped to elucidate the relationship between Navier-Stokes flow and Darcy flow, which was first derived as a phenomelogical equation.

.. tip::
    There is a great pdf on the relationship between Navier-Stokes and Darcy flow on the webpages of `Cyprien Soulaine <https://www.cypriensoulaine.com/publications-copy>`_.

During this lecture we will study single-phase hydrothermal flow in submarine hydrothermal systems. The respective solver is named :code:`HydrothermalSinglePhaseDarcyFoam`. The hydrothermal fluid flow is governed by Darcy's law (Eqn. :eq:`eq:darcy`), 
mass continuity (Eqn. :eq:`eq:conti`) and energy conservation (Eqn. :eq:`eq:temperature`) equations shown below, 

.. math::
    :label: eq:darcy 

    \vec{U} = - \frac{k}{\mu_f} (\nabla p -\rho \vec{g})
    
.. math::
    :label: eq:conti
    
    \varepsilon \frac{\partial \rho_f}{\partial t} + \nabla \cdot (\vec{U} \rho_f)

.. math::
    :label: eq:pressure
    
    \varepsilon \rho_f \left( \beta_f \frac{\partial p}{\partial t} - \alpha_f \frac{\partial T}{\partial t} \right) = \nabla \cdot \left( \rho_f \frac{k}{\mu_f} (\nabla p - \rho_f \vec{g}) \right)

.. math::
    :label: eq:temperature
    
    (\varepsilon \rho_f C_{pf} + (1-\varepsilon)\rho_r C_{pr})\frac{\partial T}{\partial t} = \nabla \cdot (\lambda_r \nabla T) - \rho_f C_{pf} \vec{U}\cdot \nabla T + \frac{\mu_f}{k} \parallel \vec{U} \parallel ^2 - \left( \frac{\partial ln \rho_f}{\partial ln T} \right)_p \frac{Dp}{Dt}

where the pressure equation :eq:`eq:pressure` is derived from continuity equation :eq:`eq:conti` and Darcy's law :eq:`eq:darcy`.

Implementation
--------------
The details of the OpenFoam implementation can be found in the |foam| documentation. Here we only show a brief summarry. :numref:`fig:htf_solution` shows how the energy equation is solved within the OpenFoam framework.


.. figure:: /_figures/solution_algorithm.*
   :align: center
   :name: fig:htf_solution

   Implementation of the energy conservation equation.

Equation-of-state
-----------------
The fluid properties like density, viscosity, specific heat are determined from the equaion-of-state of pure water. :numref:`fig:phase_diagram` shows the phase diagram of pure water. At sub-critical conditions (P< 22 MPa), the boiling curve divides the regions of liquid water and water vapor. At super-critical conditions, there is a gradual transition from a liquid-like to a vapor-like fluid phase. |foam| is a single phaes code and can only be used in regions, where a single fluid phase is present, i.e. under pure liquid water, water vapor, or supercritical conditions; boiling cannot be resolved. As we will find out later, the thermodynamic properties of water have first order control on flow dynamics and upflow temperatures in submarine hydrothermal systems. 


.. figure:: /_figures/PhaseDiagram.*
   :align: center
   :name: fig:phase_diagram

   Phase diagram of pure water.
