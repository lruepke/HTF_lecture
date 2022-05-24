.. include:: /include.rst_

.. _L04_PostProcess:

Upflow temperatures in submarine hydrothermal systems
=====================================================

During the last session we have explored free hydrothermal convection in a porous medium heated from below. That 2-D setting can be thought of as hydrothermal convection in an along axis slice above a stable melt lens - for example at a fast spreading ridge. Two things are striking about our initial results. First, in the homogenous permeability case, vent temperatures never exceed 400°C although the bottom boundary condition was set to 800°C. Shouldn't it be the hottest and most buoyant fluids that rise in free convection? Second, in the layered case, mixing at the interface to the higher permeability layer resulted in a reduction in vent temperature. What is controlling this reduction in the upwelling temperature?

If we turn to observations to adress these question, we notice another striking feature of submarine hydrothermal systems. Global compilations of vent fluid exit temperatures demonstrate that these never exceed approx. 400°C (:numref:`fig:vent_T_fig` ). 


.. figure:: /_figures/vent_T_deoth.*
   :align: center
   :name: fig:vent_T_fig
   :figwidth: 70%

   Figure taken from :cite:`Nakamura2014` showing measured vent fluid temperatures.

So, what controls the upflow temperature in submarine hydrothermal systems? This question was adressed by :cite:`jupp2000thermodynamic` in a landmark paper. They showed that the thermodynamic properties of water control the upflow temperature and can explain the observed upper limit of vent fluid temperatures. In this lecture we will have a detailed look into the underlying mechanisms.


Mechanism to limit black smoker temperature 
-----------------------------------------------

We recommend reading the original papers by :cite:`jupp2000thermodynamic` and :cite:`Jupp2004`; here we only provide a very short summary. The underlying idea is to determine from the pure water equation-of-state the optimum temperature for buoyant upwelling. Within the energy conservation equation, advective heat transport is described by a divergence term:

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

:cite:`jupp2000thermodynamic` called the term within the large brackets, :math:`F=\left(\frac{(\rho_0 - \rho) \rho h}{\mu}\right)` , fluxibility. Fluxibility is a function of fluid properties only and those properties are pressure and temperature dependent. To first order, advective heat transport is maximized, where :math:`\nabla \cdot F` is maximum and this divergence can be approximated within the bottom thermal boundary as :math:`\nabla \cdot F \approx \frac{\partial}{\partial z} F \approx \frac{\partial}{\partial T} F` . 


:numref:`fig:Fluxibility_Water` a shows fluxibility as a function of pressure and temperature. It is clear that this function has a distinct peak at temperature of approximately 400°C. :numref:`fig:Fluxibility_Water` b shows sections of constant pressure and also the derivative of F with temperature. The peaks in these functions mark the temperature for which buoyant heat transport is most efficient for a given pressure. Later we will explore this further in numerical convection experiments.

.. figure:: /_figures/Fluxibility_Water.*
   :align: center
   :name: fig:Fluxibility_Water

   Fluxibility F as a function of temperature and pressure (a), profile of F along temperature with constant pressure (b).


Let's have a look into the controlling fluid properties density, specific enthalpy, and dynamic viscosity in P,T space.

.. only:: html

    .. tab:: P-T H2O phase diagram

        .. figure:: /_figures/sphx_glr_plot_PhaseDiagram_pT_003_2_0x.*
            :align: center

            Phase regions of pure water in P-T space. We are mainly interested in the supercritical region above the critical end point, when water is always in the single phase regime.

    .. tab:: Density

        .. figure:: /_figures/sphx_glr_plot_Density_pT_001_2_0x.*
            :align: center
            :figwidth: 60%

            Density (:math:`Kg/m^3`) as a function of pressure and temperature. There is a monotonic rise in density with increasing pressure and and decreasing temperature - just as intuition would tell us. 

    .. tab:: Specific enthalpy

        .. figure:: /_figures/sphx_glr_plot_Enthalpy_pT_001_2_0x.*
            :align: center
            :figwidth: 60%

            Specific enthalpy (:math:`J/Kg`) as a function of pressure and temperature. The energy "stored" in a kilogram of water rises with increasing tempeature - also no surprises here. 

    .. tab:: Dynamic viscosity

        .. figure:: /_figures/sphx_glr_plot_Viscosity_pT_001_2_0x.*
            :align: center
            :figwidth: 60%

            Viscosity (:math:`Pa s`) as a function of pressure and temperature. Here things get interesting; the viscosity is first going down with increasing temperature (as expected) but then increases again within the vapor phase field. 


Looks like viscosity may play a big role in determining flow dynamics and upwelling temperatures. Let's explore this in more detail.