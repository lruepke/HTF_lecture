.. include:: /include.rst_

.. _L05_Exercise2:

Exercise 2: Fluxibility and Tvent
=================================

Before going into more detail, let's have a closer into the case we explored last week. 

    * Load the homogeneous or layered case with constant
    * Use the paraview calculator to create a new variable T_Celsius (by doing T_Celsius = T - 273.15)
    * Make a z-normal slice
    * And explore the 100,200,300,400,500 °C contours through time?

.. figure:: /_figures/Fluxibility_postprocess.*
   :align: center
   :name: fig:Fluxibility_postprocess_fig

   Notice how the diffusive thermal boundary layer becomes unstable between the 300°C and 400°C isotherms. This is what concept of fluxibility tells us!


:cite:`jupp2000thermodynamic` concept of fluxibility was later extended to explore interrelations between permeability on vent temperatures :cite:`driesner2010`, the structure of hydrothermal flow zone :cite:`coumou2008`, and fault-controlled off-axis hydrothermal system :cite:`andersen2015`. Some of these papers are available in the e-learning site of this course for registered CAU students.

Let's look into the original setup! We will  try to reproduce the results presented by :cite:`jupp2000thermodynamic`. The setup is very similar to the simple single plume driven by a gaussian heat source we had simulated previously. 

Model setup and running the case
--------------------------------

We are simulating a 3700 x 1000m sized box with a gaussian constant temperature boundary condition at the bottom. The general model setup is shown in the Figure below. The OpenFoam case can be downloaded from here: :download:`Jupp_Schultz <cases/Jupp_Schultz.zip>`.

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

Download the case file, unzip it and copy it into your shared folder that is accessible both from your local system and from the docker container. Check out the usual files like :code:`system/blockMeshDict`, :code:`run.sh`, :code:`0/T`. Make sure you understand the setup!

Code modifications and local Rayleigh number
--------------------------------------------

We have learned about "fluxibility" but let's go one step back and first think about the vigor of convection in general. The tendency of a fluid to convect and the vigor/regime of free convection is often expressed by a dimensionless number, the famous Rayleigh number, :math:`Ra`. It naturally shows up as the driving force of convection, when the governing equations are non-dimensionalized, or one can follow a more classic derivation. There is a good explanation on `wikipedia <https://en.wikipedia.org/wiki/Rayleigh_number>`_ . Formally :math:`Ra` is the ratio between the time scales of heat transport by diffusion and convection. 

Let's have a quick look and do a quick dimensional analysis. Simple heat diffusion looks like this:

.. math::
    :label: eq:dim_ana_1

    \frac{\partial T}{\partial t} = \kappa \frac{\partial^2 T}{\partial x^2}


Let's make the variables dimensionless by introducing characteristic scales:

.. math::
    :label: eq:dim_ana_2

    T' = \frac{T}{\Delta T}

.. math::
    :label: eq:dim_ana_3

    x' = \frac{x}{L}

.. math::
    :label: eq:dim_ana_4

    t' = \frac{t}{\frac{L^2}{\kappa}}

There we have it. The characteristic time of diffusion is :math:`\tau_{diff} =\frac{L^2}{\kappa}` with units of :math:`s`.

Alright, let's look at convection. In Darcy flow the velocity is proportional to the density difference according to :math:`u \sim \frac{\Delta \rho g k}{\mu}`. If we express the density difference via the thermal expansion and divide the characteristic length, :math:`L` by this velocity, we get :math:`\tau_{conv} =\frac{L \mu}{\Delta T \alpha \rho g k}`, which again has units of :math:`s`.

So the Rayleigh number is:

.. math::
    :label: eq:dim_ana_5

    Ra = \frac{\frac{L^2}{\kappa}}{\frac{L \mu}{\Delta T \alpha \rho g k}}

.. math::
    :label: eq:dim_ana_6
    
    Ra = \lvert \frac{\rho \Delta T \alpha g k L}{\kappa \mu} \rvert


With this analysis in mind, we can check local Rayleigh numbers and what we expect is that right where the fluxibity is highest, we should also get high Rayleigh numbers indicating that the fluid wants to start convecting. The Rayleigh number derived above is a measure for the vigor of convection in the entire domain. We can follow :cite:`jupp2000thermodynamic` and derive a local Rayleigh number by relating convective and diffusive fluxes:

.. math::
    :label: eq:dim_ana_6
    
    Ra_L = \lvert \frac{\nabla \cdot (\vec{u} \rho h)}{\lambda \nabla^2 T} \rvert


The convective energy transport, :math:`\nabla \cdot (\vec{u} \rho h)` term is back! If you read :cite:`jupp2000thermodynamic` carefully, you find a dimensional analysis that shows that :math:`Ra_L` becomes maximum where that transport term is maximum and we are back to the fluxibility concept.


Enough theroy, let's put this into a model. The compute local Rayleigh numbers, we have to add something to the code. There is a dynamic code section within :code:`system/conrolDict` that allows outputting more variables including the thermodynamic properties of water. Look how easy it is to evaluate the local convective and diffusive fluxes.

.. code-block:: foam
    :linenos:

    functions
    {
        calPlumeT
        {
            libs                ("libutilityFunctionObjects.so");
            type                coded;
            enabled             true;
            writeControl        adjustableRunTime;
            writeInterval       $writeInterval;
            name                plumeTemperature;
            codeWrite
            #{
                //get maximum tempeature on the top boundary
                label patchID = mesh().boundaryMesh().findPatchID("top"); 
                const volScalarField& T = mesh().lookupObject<volScalarField>("T");
                const volScalarField& mu = mesh().lookupObject<volScalarField>("mu");
                const volScalarField& Cp = mesh().lookupObject<volScalarField>("Cp");
                const volScalarField& rho = mesh().lookupObject<volScalarField>("rho");
                const surfaceScalarField& phi = mesh().lookupObject<surfaceScalarField>("phi");
                const volScalarField& h = mesh().lookupObject<volScalarField>("enthalpy");
                double kr = 2.0;
                // calculate local Rayleigh number 
                volScalarField Ra_L("LocalRayleigh",mag(fvc::div(phi, h) / (kr*fvc::laplacian(T))));
                // write properties to output files
                mu.write(); Cp.write(); h.write(); Ra_L.write();
                // write vent temperature
                std::ofstream fout("ventT.txt",std::ofstream::app);
                // Info<<"Plume Temperature: "<<mesh().time().value()/31536000<<"\t"<<Foam::gMax(T.boundaryField()[patchID])-273.15<<" C"<<endl;
                fout<<mesh().time().value()/31536000<<"\t"<<Foam::gMax(T.boundaryField()[patchID])-273.15<<std::endl;
                fout.close();
            #};
        }
    }


That code block requires the specific enthalpy to be accessible by the code. The version of |foam| we are using does not have that activated. Therefore we need to update it to the newest version. Thankfully that's easy to do; go into the docker container shell and

.. code-block:: bash

    cd $HOME
    ./getHydrothermalFoam_latest.sh

This will update the source code of |foam| and recompile it. Just wait for it to complete. Now you can run the case!

It's possible that we have been a bit overambitious in our chosen numerical resolution. In case the run-time is too long, just reduce the horizontal and vertical resolution. Make the mesh with :code:`blockMesh` and call the solver; or execute the :code:`run.sh` script. If you are feeling courages today, you can also try our running it in parallel! Check out the :code:`run_par.sh` script for that.
You can also shorten the total run-time (to say 40 yrs) in :code:`system/controlDict`; we are mainly interested in the early stages anyway. 


Post-processing
---------------

Plotting temperature
^^^^^^^^^^^^^^^^^^^^

After running the case, check it out in paraview and make sure you understand the results and that everything looks fine. While paraview is great for 3-D processing, sometimes python can be more powerfull in detailed post-processing. 

The figure below shows an example plot of the final temperature solution. It basically does the same thing as paraview, but let's walk through it. 

.. figure:: /_figures/T_Jupp_Schultz_steady.*
   :align: center
   :name: fig:T_Jupp_Schultz_steady

   The steady-state temperature distribution. See :ref:`/lectures/L05/cases/Jupp_Schultz/jupyter/Plot_CaseResults.ipynb#Temperature-+-velocity-field-(quiver)` in the notebook.


Click on the link in the figure caption and go through the notebook. Safest way of getting this to work is to copy the steps into local notebooks (see :ref:`Installation guide`). 


You can also start up paraview and explore everything there! Check which temperature contours become unstable first, plot vectors and local Rayleigh numbers.