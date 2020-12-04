.. include:: /include.rst_

.. _L04_Exercise2:

Excercise 2
===========

We will now try to reproduce the results presented by :cite:`jupp2000thermodynamic`, the case files :download:`Jupp_Schultz <cases/Jupp_Schultz.zip>`.

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


.. figure:: /_figures/T_Jupp_Schultz_steady.*
   :align: center
   :name: fig:T_Jupp_Schultz_steady

   The steady-state temperature distribution. See :ref:`/lectures/L04/cases/Jupp_Schultz/jupyter/Plot_CaseResults.ipynb#Temperature-+-velocity-field-(quiver)` in the notebook.

.. figure:: /_figures/Jupp_Schultz_early.*
   :align: center
   :name: fig:Jupp_Schultz_early

   Early stages of the simulation. See :ref:`/lectures/L04/cases/Jupp_Schultz/jupyter/Plot_CaseResults.ipynb#Early-stage` in the notebook.
   



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

.. code-block:: foam
    :linenos:
    :name: lst:controlDict:functions:props_Ral
    :caption: Post-processing function in controlDict to output properties and calculate local Rayleigh number.

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

.. toctree::
    :maxdepth: 2

    cases/Jupp_Schultz/jupyter/Plot_CaseResults.ipynb