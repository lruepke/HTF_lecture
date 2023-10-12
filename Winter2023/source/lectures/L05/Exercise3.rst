.. include:: /include.rst_

.. _L05_Exercise3:

Exercise 3: Permeability and Tvent
==================================

We now understand a bit better why the approx. 400°C fluids are rising towards the surface and the hotter fluids are remaining stagnant within the reaction zone. Let's do the next step and understand what is happening at the interface between lower and higher permeabilities. 

Go back to our two-layered case and post-process the solution in paraview a bit more. 

    * Load the layered case with constant
    * Use the paraview calculator to create a two new variables T_Celsius (by doing T_Celsius = T - 273.15) and mass flux (U*rho)
    * Make a z-normal slice
    * Puzzle temperature contours and vertical mass flux together and explore what happens at the interface.

.. figure:: /_figures/Flux_interface.*
   :align: center
   :name: fig:Flux_interface_fig


Clearly, mass flux is increasing and temperature is decreasing. What do you think is happening to the convective heat flux?

We will now explore the interplay between permeability and bottom heat flux. The research question is: what kind of permeability is needed to sustain a high-temperature hydrothermal system? This question has been addressed by :cite:`driesner2010interplay`.

Let's first go through the theory. The flux into the reaction zone of hydrothermal convection cell can be written as:

.. math::
    :label: eq:heat_flux_rz

    \vec{q_{in}} = 2\frac{L\lambda (T_D - T_U)}{H_R}

:math:`L` is the half-width of the reaction zone, :math:`\lambda` is the thermal conductivity, :math:`T_D` the temperature of the driving heat source, :math:`T_U` the upflow temperature, and :math:`H_R` the thickness of the reaction zone.

Let's follow :cite:`jupp2000thermodynamic` and explore the heat transported by hydrothermal convection. For this we have to make a few assumptions. The first is that the pressure gradient in the upwflow is close to cold hydrostatic :math:`\rho_0 g`, which is a reasonable assumption based on convection experiments. 

Now we can use Darcy's law to spell out the vertical velocity and mass flux:

.. math::
    :label: eq:mass_flux_1

    \vec{u} = -\frac{k}{\mu}( \nabla P - \rho \vec{g} )

.. math::
    :label: eq:mass_flux_2

    \frac{\partial P}{\partial z} = \rho_0 g_z

.. math::
    :label: eq:mass_flux_3

    u_z = -\frac{k}{\mu_U}g_z(\rho_0 - \rho_U)

.. math::
    :label: eq:mass_flux_4

    \rho_U u_z = -\frac{k}{\mu_U} \rho_U g_z(\rho_0 - \rho_U)

with :math:`_U` always referring to properties within the upflow zone, which has a temperature of :math:`T_U`.

We can spell out the heat output by multiplying this with the difference in specific enthalpy :math:`h [J/Kg]` between the upflow (:math:`h_U`) and the recharge zone (:math:`h_0`). To get the total heat flux, we also need to multiply with two times the half-width (:math:`L`) of the upflow zone:

.. math::
    :label: eq:mass_flux_5

    2L\rho_U u_z (h_U - h_0) = 2gk \left[ \frac{\rho_U (h_U - h_0) (\rho_0 - \rho_U) }{\mu_U}\right]L


And by setting this equal to the heat input, we get:

.. math::
    :label: eq:heat_balance

    2\frac{L\lambda (T_D - T_U)}{H_R} = 2gk \left[ \frac{\rho_U (h_U - h_0) (\rho_0 - \rho_U) }{\mu_U}\right]L


Note how we "flipped" the sign of g for better readability and omitted the :math:`T_U` index. If we replace the right-hand side with a given heat flux :math:`Q` we get:

.. math::
    :label: eq:heat_balance_2

    Q = 2gk \left[ \frac{\rho_U (h_U - h_0) (\rho_0 - \rho_U) }{\mu_U}\right]L


Nice, we have an equation that relates heat input, permeability, and upflow temperature. Unfortunately, it is not that useful as the upflow temperature is hidden as an implicit term in the fluid properties. 

Try it out!
-----------

Let's implement the equation in python and explore the solution and implications. Go through the attached notebook, fill out the missing pieces and think about what this means for upflow temperatures in genetral and for heterogeneous rock in particular.

.. toctree::
    :maxdepth: 2

    driesner_2010.ipynb


Numerical solution
-------------------

Next we will explore this using our numerical model.  A basic setup can be downloaded from here: :download:`Driesner2010 <cases/Driesner2010.zip>`.

The magma heat source is simulated by a heat flux boundary condition, which can be set by customized boundary condition type of :code:`hydrothermalHeatFlux` (see `doc <https://www.hydrothermalfoam.info/manual/en/Models_Equations/index.html#heat-flux-bc-hydrothermalheatflux>`_ ).
The model in :cite:`driesner2010interplay` is 1 m thick, 3 km wide, 1 km height, the heat source is simulated by a Gaussian-shaped heat flux profile with total heat input 86 kW/m, half-width 500 m and center 1500 m. 

Let's see how to set this boundary condition in HydrothermalFoam. 
:code:`hydrothermalHeatFlux` supports Gaussian-shape (:code:`shape gaussian2d;`) distribution and the shape is defined by four parameters :code:`x0, qmax, qmin, c`, the equation of the shape is 

.. math::
    :label: eq:gauss_hf

    q_h(x) = q_{min} + (q_{max}-q_{min})e^{-\frac{(x-x_0)^2}{2c^2}}

for a normal Gaussian-shape profile in this model, :code:`qmin` is set to zero and :code:`x0` is set to 1500. 
From equation :eq:`eq:gauss_hf`, we can get the half-width and total heat flux (approximately) are :math:`c\sqrt{2ln2}` and :math:`(q_{max}-q_{min})\sqrt{2\pi}c`, respectively. According to the model setup, it's easy to get :math:`c = \frac{500}{\sqrt{2ln2}} = 424.661` and :math:`q_{max} = \frac{86}{c\sqrt{2\pi}} = 0.0808\ kW/m^2`.
There the temperature boundary condition at the bottom patch can be set as,

.. tab:: Code snippet

    .. code-block:: foam
        :caption: Temperature boundary condition of bottom patch
        :name: lst:BC_hydrothermalHeatFlux

        bottom
        {
            type        hydrothermalHeatFlux;
            q           uniform 0.05; //placeholder
            value       uniform 0; //placeholder
            shape       gaussian2d;
            x0          1500;
            qmax        80.8;
            qmin        0;
            c           424.661;
        }

.. tab:: Heat flux profile 

    .. plot:: 

        xmin,xmax=0,3000
        width_half,q_total = 500, 86000
        c=width_half/np.sqrt(2*np.log(2))
        qmax,qmin,x0=q_total/(c*np.sqrt(2*np.pi)),0,1500
        fig=plt.figure(figsize=(6,3))
        ax=plt.gca() 
        x=np.linspace(xmin,xmax,1000)
        q=qmin+(qmax-qmin)*np.exp(-(x-x0)*(x-x0)/(2*c*c))
        ax.plot(x,q)
        ax.fill_between(x,q,0,color='lightgreen',alpha=0.5)
        ax.hlines(0.5, 0,0.5,transform=ax.transAxes, ls='dashed',lw=0.5, color='gray')
        ax.vlines(0.5, 0.5,1,transform=ax.transAxes, ls='dashed',lw=0.5, color='gray')
        ax.plot([x0-width_half,x0],[qmax/2, qmax/2],color='r',lw=2)
        ax.text(x0-width_half/2, qmax/2, 'half width\n%.0f'%(width_half),va='top',ha='center')
        q_total_cal=0
        dx=(x.max()-x.min())/(len(x)-1)
        # integrate heat flux along profile
        for i in range(0,len(x)):
            q_total_cal = q_total_cal + dx*q[i]
        ax.text(x0,qmin*1.2,'Integrated total heat flux\n%.1f kw'%(q_total_cal/1000),va='bottom',ha='center',fontweight='bold')
        text_parms='$x_{0}$ = %.0f\n$q_{min}$ = %.0f\nc = %.1f m\n$q_{max}$ = %.1f w/m$^2$'%(x0,qmin,c,qmax)
        ax.text(0.98,0.98,text_parms,va='top',ha='right',transform=ax.transAxes)
        ax.set_ylabel('Heat flux (W/m$^{\mathregular{2}}$)')
        ax.set_xlabel('Distance (m)')
        ax.set_title('Gaussian shape heat flux profile')
        ax.set_ylim(q.min(),q.max())
        plt.tight_layout()
        plt.show()


Excercise
---------

Set up a new run, for example with :math:`k=5e^{-14} m^2`, and compare the vent temperatures of the numerical simulation with the predictions of our semi-analytical solution.