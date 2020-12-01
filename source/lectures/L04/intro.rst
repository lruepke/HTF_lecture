.. include:: /include.rst_

Porous Flow - Submarine Hydrothermal Systems
============================================


.. _Jupp_Schultz:

Mechanism to limit black smoker temperature 
-----------------------------------------------

Download the case files :download:`Jupp_Schultz <cases/Jupp_Schultz.zip>`



Properties calculation details can be found in the
 
.. figure:: /_figures/Properties_Water.*
   :align: center
   :name: fig:Properties_Water

   Water properties as a function of temperature and pressure, using python package of iapws_.

.. figure:: /_figures/Fluxibility_Water.*
   :align: center
   :name: fig:Fluxibility_Water

   Fluxibility F as a function of temperature and pressure (a), profile of F along temperature with constant pressure (b).

.. plot::

    xmin,xmax=-1700,1700
    Tmin,Tmax=20,1000
    wGauss,x0=10,0
    fig=plt.figure(figsize=(6,3))
    ax=plt.gca()
    x=np.linspace(xmin,xmax,500)
    y=Tmin+(Tmax-Tmin)*np.exp(-(x-x0)*(x-x0)/(2*wGauss*wGauss))
    ax.plot(x,y)
    ax.set_ylabel('Temperature ($^{\circ}$C)')
    ax.set_xlabel('Distance (m)')
    ax.set_title('Gaussian shape temperature boundary condition')
    ax.set_ylim(y.min(),y.max())
    plt.tight_layout()
    plt.show()