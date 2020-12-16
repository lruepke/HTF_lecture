.. include:: /include.rst_

.. _L04_Exercise2:

Excercise 2
===========

:cite:`jupp2000thermodynamic` published a landmark paper on how A landmark paper on how the thermodynamic properties of water determine temperatures in hydrothermal upflow zones. Their concept of "fluxibility" was later extended to explore interrelations between permeability on vent temperatures :cite:`driesner2010`, the structure of hydrothermal flow zone :cite:`coumou2008`, and fault-controlled off-axis hydrothermal system :cite:`andersen2015`. Some of these papers are available in the e-learning site of this course for registered CAU students.

.. admonition:: Objectives of this exercise

    - understand how thermodynamic properties control upflow temperatures
    - learn how to use python to post-process OpenFoam cases

We will now try to reproduce the results presented by :cite:`jupp2000thermodynamic`. The setup is very similar to the simple single plume driven by a gaussian heat source we had simulated previously. 

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

There is a dynamic code section within :code:`system/conrolDict` that allows outputting more variables including the thermodynamic properties of water.

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


Post-processing
---------------

Plotting temperature
^^^^^^^^^^^^^^^^^^^^

After running the case, check it out in paraview and make sure you understand the results and that everything looks fine. While paraview is great for 3-D processing, sometimes python can be more powerfull in detailed post-processing. 


.. tip::
    Have you noticed the python terminal in paraview? Combining python with paraview can be extremely powerful! If you are curious, check it out on the internet. We will also show some examples later in the course.


The figure below shows an example plot of the final temperature solution. It basically does the same thing as paraview, but let's walk through it. 

.. figure:: /_figures/T_Jupp_Schultz_steady.*
   :align: center
   :name: fig:T_Jupp_Schultz_steady

   The steady-state temperature distribution. See :ref:`/lectures/L04/cases/Jupp_Schultz/jupyter/Plot_CaseResults.ipynb#Temperature-+-velocity-field-(quiver)` in the notebook.


Click on the link in the figure caption and go through the notebook. Safest way of getting this to work is to copy the steps into local notebooks (see :ref:`Installation guide`). 

.. code-block::

    code plot_temp.ipynb

Copy the sections into your local notebook, try to execute the sections (press SHIFT+RETURN), add some markdown sections with comments. 
   

Exploring the parameter space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Effect of cell size (**CS**) and full width at half maximum (**FWHM**),

.. tab:: CS=10,FWHM=300

    .. figure:: /_figures/T_Jupp_Schultz300_steady.*
        :align: center

        Cell size = 10 m, FWHM = 300 m.

.. tab:: CS=10,FWHM=400

    .. figure:: /_figures/T_Jupp_Schultz400_steady.*
        :align: center

        Cell size = 10 m, FWHM = 400 m.

.. tab:: CS=2.5,FWHM=300

    .. figure:: /_figures/T_Jupp_Schultz300_2.5_steady.*
        :align: center

        Cell size = 2.5 m, FWHM = 300 m.


.. figure:: /_figures/Jupp_Schultz_early.*
   :align: center
   :name: fig:Jupp_Schultz_early

   Early stages of the simulation. See :ref:`/lectures/L04/cases/Jupp_Schultz/jupyter/Plot_CaseResults.ipynb#Early-stage` in the notebook.


.. toctree::
    :maxdepth: 2

    cases/Jupp_Schultz/jupyter/Plot_CaseResults.ipynb