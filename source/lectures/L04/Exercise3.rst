.. include:: /include.rst_

.. _L04_Exercise3:

Excercise 3
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

    q_h(x) = q_{min} + (q_{max}-q_{min})e^{-\frac{(x-x_0)^2 + (y-z_0)^2}{2c^2}}

for a normal Gaussian-shape profile in this model, :code:`qmin` is set to zero and :code:`x0` is set to 1500. 
From equation :eq:`eq:gauss_hf`, we can get the half-width and total heat flux (approximately) are :math:`c\sqrt{2ln2}` and :math:`(q_{max}-q_{min})\sqrt{2\pi}c`, respectively. According to the model setup, it's easy to get :math:`c = \frac{500}{\sqrt{2ln2}} = 424.661` and :math:`q_{max} = \frac{86}{c\sqrt{2\pi}} = 0.0808\ kw/m^2`.
There the temperature boundary condition at the bottom patch can be set as,

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


.. figure:: /_figures/ventT_time_perm.*
   :align: center
   :name: fig:ventT_time_perm

   :ref:`/lectures/L04/cases/Driesner2010/jupyter/Plot_CaseResults.ipynb#Maximum-vent-temperature-changes-with-time`.

.. figure:: /_figures/ventT_perm.*
   :align: center
   :name: fig:ventT_perm
   
   :ref:`/lectures/L04/cases/Driesner2010/jupyter/Plot_CaseResults.ipynb#Vent-temperature-as-a-function-of-permeability`.

.. toctree::
    :maxdepth: 2

    cases/Driesner2010/jupyter/Plot_CaseResults.ipynb