.. include:: /include.rst_

.. _L04_Exercise3:

Exercise 3
===========

.. warning::
   This section is under construction!

We will now try to reproduce the results presented by :cite:`driesner2010interplay`, the case files :download:`Driesner2010 <cases/Driesner2010.zip>`.

The magma heat source is simulated by a heat flux boundary condition, which can be set by customized boundary condition type of :code:`hydrothermalHeatFlux` (see `doc <https://www.hydrothermalfoam.info/manual/en/Models_Equations/index.html#heat-flux-bc-hydrothermalheatflux>`_ ).
The model in :cite:`driesner2010interplay` is 1 m thick, 3 km wide, 1 km height, the heat source is simulated by a Gaussian-shaped heat flux profile with total heat input 86 km, half-width 500 m and center 1500 m. 

Let's see how to set this boundary condition in HydrothermalFoam. 
:code:`hydrothermalHeatFlux` supports Gaussian-shape (:code:`shape gaussian2d;`) distribution and the shape is defined by four parameters :code:`x0, qmax, qmin, c`, the equation of the shape is 

.. math::
    :label: eq:gauss_hf

    q_h(x) = q_{min} + (q_{max}-q_{min})e^{-\frac{(x-x_0)^2}{2c^2}}

for a normal Gaussian-shape profile in this model, :code:`qmin` is set to zero and :code:`x0` is set to 1500. 
From equation :eq:`eq:gauss_hf`, we can get the half-width and total heat flux (approximately) are :math:`c\sqrt{2ln2}` and :math:`(q_{max}-q_{min})\sqrt{2\pi}c`, respectively. According to the model setup, it's easy to get :math:`c = \frac{500}{\sqrt{2ln2}} = 424.661` and :math:`q_{max} = \frac{86}{c\sqrt{2\pi}} = 0.0808\ kw/m^2`.
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

.. figure:: /_figures/ventT_time_perm.*
   :align: center
   :name: fig:ventT_time_perm

   :ref:`/lectures/L05/cases/Driesner2010/jupyter/Plot_CaseResults.ipynb#Maximum-vent-temperature-changes-with-time`.

.. figure:: /_figures/ventT_perm.*
   :align: center
   :name: fig:ventT_perm
   
   :ref:`/lectures/L05/cases/Driesner2010/jupyter/Plot_CaseResults.ipynb#Vent-temperature-as-a-function-of-permeability`.

.. toctree::
    :maxdepth: 2

    cases/Driesner2010/jupyter/Plot_CaseResults.ipynb