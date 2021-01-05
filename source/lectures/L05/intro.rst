.. include:: /include.rst_

Structural controls on hydrothermal flow
=====================================================

A striking observations from global vent fluid measurements is that fluid exit temperatures apparently never exceed 400°C (:numref:`fig:vent_T_fig` ). 


.. figure:: /_figures/vent_T_deoth.*
   :align: center
   :name: fig:vent_T_fig
   :figwidth: 70%

   Figure taken from :cite:`Nakamura2014` showing measured vent fluid temperatures.

Maybe you have already noticed in our convection experiments in the previous session that the upflow temperatures are always much lower than the bottom boundary condition implying that the hottest fluids are actually stagnant and do not rise towards the seafloor. :cite:`jupp2000thermodynamic` published a landmark paper in which they showed that the thermodynamic properties of water control the upflow temperature and can explain the observed upper limit of vent fluid temperatures. In this lecture we will have a detailed look into the proposed mechanisms.


Mechanism to limit black smoker temperature 
-----------------------------------------------

We recommend reading the original papers by :cite:`jupp2000thermodynamic` and :cite:`Jupp2004`; here we only provide a very short summarry. The underlying idea is to determine from the pure water equation-of-state the optimum temperature for buoyant upwelling. Within the energy conservation equation, advective heat transport is described by a divergence term:

.. math::
    :label: eq:div_e
    
    \nabla \cdot (\rho h \vec{U})

Following :cite:`jupp2000thermodynamic`, we can assume that the pressure gradient, driving upflow, can be approximated as cold hydrostatic so that we can express the upflow velocity as:

.. math::
    :label: eq:upflow_js
    
    U_z = \frac{k}{\mu} \left( \nabla P - \rho g_z \right) \approx \frac{k}{\mu}(\rho_0 g_z - \rho g_z) = g_z \frac{k}{\mu} (\rho_0 - \rho)

Now we plug this expression into :eq:`eq:div_e` to get:

.. math::
    :label: eq:udiv_e2
    
    \nabla \cdot (\rho h \vec{U}) = \frac{\partial}{\partial z} \left(\frac{(\rho_0 - \rho) \rho h}{\mu}\right) g_z k

:cite:`jupp2000thermodynamic` called the term wihtin the large brackets, :math:`F=\left(\frac{(\rho_0 - \rho) \rho h}{\mu}\right)` , fluxibility. Fluxibility is a function of fluid properties only and those properties are pressure and temperature dependent. To first order, advective heat transport is maximized, where :math:`\nabla \cdot F` is maximum and this divergence can be approximated within the basal boundary as :math:`\nabla \cdot F \approx \frac{\partial}{\partial z} F \approx \frac{\partial}{\partial T} F` . 


:numref:`fig:Fluxibility_Water` a shows fluxibility as a function of pressure and temperature. It is clear that this function has a distinct peak at temperature of approximately 400°C. :numref:`fig:Fluxibility_Water` b shows sections of constant pressure and also the derivative of F with temperature. The peaks in these functions mark the temperature for which buoyant heat transport is most efficient for a given pressure. Later we will explore this further in numerical convection experiments.

.. figure:: /_figures/Fluxibility_Water.*
   :align: center
   :name: fig:Fluxibility_Water

   Fluxibility F as a function of temperature and pressure (a), profile of F along temperature with constant pressure (b).
