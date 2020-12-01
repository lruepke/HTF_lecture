.. include:: /include.rst_

Porous Flow - Submarine Hydrothermal Systems
============================================


.. _Jupp_Schultz:

Mechanism to limit black smoker temperature 
-----------------------------------------------


Download the case files :download:`Jupp_Schultz <cases/Jupp_Schultz.zip>`


.. figure:: /_figures/Properties_Water.*
   :align: center
   :name: fig:Properties_Water

   Water properties as a function of temperature and pressure, using python package of iapws_. Properties calculation details can be found in the `Jupyter-notebook <cases/Jupp_Schultz/jupyter/Plot_Lectures.ipynb>`_
 
.. figure:: /_figures/Fluxibility_Water.*
   :align: center
   :name: fig:Fluxibility_Water

   Fluxibility F as a function of temperature and pressure (a), profile of F along temperature with constant pressure (b).

The mode setup is based on :cite:`jupp2000thermodynamic`, the heat source is represented by a Gaussian shape constant temperature profile shown below.

.. plot::

    xmin,xmax=-1700,1700
    Tmin,Tmax=20,1000
    wGauss,x0=400,0
    xigma = wGauss/2.355
    # see https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    fig=plt.figure(figsize=(6,3))
    ax=plt.gca()
    x=np.linspace(xmin,xmax,500)
    y=Tmin+(Tmax-Tmin)*np.exp(-(x-x0)*(x-x0)/(2*xigma*xigma))
    ax.plot(x,y)
    ax.axvspan(x0-wGauss/2,x0+wGauss/2,ymax=0.5, alpha=0.5, color='lightgreen')
    ax.text(0.5,0.25,'FWHM\n%.0f m'%(wGauss), va='center',ha='center',transform=ax.transAxes)
    ax.set_ylabel('Temperature ($^{\circ}$C)')
    ax.set_xlabel('Distance (m)')
    ax.set_title('Gaussian shape temperature boundary condition')
    ax.set_ylim(y.min(),y.max())
    plt.tight_layout()
    plt.show()