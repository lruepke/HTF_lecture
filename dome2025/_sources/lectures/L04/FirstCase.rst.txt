.. include:: /include.rst_

.. _L04_FirstCase:

Hydrothermal convection test case
====================================

Prepare case files
------------------
To get started we will run the **Regular2DBox** case from the cookbook directory of |foam|. This cookbooks describes how we can simulate a simple hydrothermal convection cell. It resolves hydrothermal convection driven by a gaussian-shaped constant temperature boundary condition at the bottom. 

Copy the  case into your shared working directory (probably $HOME/HydrothermalFoam_runs). You need to do this within the docker container (your right-hand shell in Visual Studio Code if you followed the recommended setup). Cd into your shared folder and type this:

.. code-block:: bash 
    :name: lst:cp2dBoxToWorkDir

    cd $HOME/HydrothermalFoam_runs
    cp -r $HOME/HydrothermalFoam/cookbooks/2d/Regular2DBox . 

Check out the directory structure shown in :numref:`lst:2dbox:tree`.

.. code-block:: bash 
    :linenos:
    :emphasize-lines: 3-5
    :name: lst:2dbox:tree
    :caption: File tree structure of the Regular2DBox case.

    .
    |-- 0
    |   |-- T
    |   |-- p
    |    -- permeability
    |-- a.foam
    |-- clean.sh
    |-- constant
    |   |-- g
    |    -- thermophysicalProperties
    |-- run.sh
     -- system
        |-- blockMeshDict
        |-- controlDict
        |-- fvSchemes
         -- fvSolution.
 
The 0 directory now has entries for T (temperature) and p (pressure) our new primary variables, and for permeability, which we will discuss later.

.. tip::
    Most OpenFoam cases include scripts like :code:`run.sh` and :code:`clean.sh`. The :code:`run.sh` script is a good starting point for "understanding" a case. It lists all commands that have to be executed (e.g. meshing, setting of properties, etc.) to run a case. The :code:`clean.sh` script cleans up the case and deletes e.g. the mesh and all output directories. Have a look into these files and see if you understand them!

The 0 directory contains all initial and boundary conditions, the system folder contains all controlling parameter files, and the constant folder contains constant properties like the mesh - which we will create next.

Equation of state and thermophysical properties
------------------------------------------------
To compute the thermodynamic properties of water, we use Xthermo, a a novel implementation of the H2O-NaCl equation-of-state by :cite:`Driesner2007`. 

xThermoProperties
^^^^^^^^^^^^^^^^^^
Just like everthing else, the paramters of the equation-of-state are set in a dictionary. Within the :code:`constant` folder, there is a dictionary called :code:`xThermoProperties`. It's structure is shown in :numref:`lst:2dbox:xThermoProps`.


.. code-block:: foam 
    :name: lst:2dbox:xThermoProps
    :linenos:
    :emphasize-lines: 17, 18, 22


    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  5                                     |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      xThermoProperties;
    }

    fluid H2O-NaCl; //H2O, H2O-NaCl
    backend IAPS84; //IAPS84, IAPWS95, IAPWS95_CoolProp

    H2O-NaCl
    {
        constX 0.0;
    }

The key properties are the choice of fluid (H2O or H2O-NaCl) and the backend (IAPS84, IAPWS95, IAPWS95_CoolProp). The backend defines the implementation of the equation-of-state. The constX parameter sets the constant salinity of the fluid in mass fraction (no phase separation phenomena are implemented in this version, so that the salinity is a constant). Here we set it to zero, so that we simulate pure water.


transportProperties
^^^^^^^^^^^^^^^^^^^^
Next we can review the properties of the solid matrix, which are set within the  :code:`constant/transportProperties` dictionary. Its structure is shown in :numref:`lst:2dbox:transProps`.

.. code-block:: foam
    :name: lst:2dbox:transProps

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  5                                     |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      transportProperties;
    }

    porosity porosity [0 0 0 0 0 0 0] 0.1;
    kr kr [1 1 -3 -1 0 0 0] 1.5;
    cp_rock cp_rock [0 2 -2 -1 0 0 0] 880;
    rho_rock rho_rock [1 -3 0 0 0 0 0] 3000;



Mesh generation
---------------

The case is run on a simple 2-d-box-like geometry and the mesh is build using :code:`blockMesh`, just like in the previous lecture on cavity flow. Look at :code:`blockMeshDict` and check that you sill understand the structure. Afterwards, you can create the mesh:

.. code-block:: bash 

    blockMesh

After making the mesh, you can use Paraview_ to visualize it,

.. code-block:: bash

    touch a.foam
    paraview a.foam 

        
Boundary conditions
-------------------

Next we need to set boundary conditions. Open the file T inside the 0 directory from your local left-hand shell.

.. code-block:: bash 

    code 0/T

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 17, 29,35,41-52
    :name: lst:2dbox:bc
    :caption: Boundary conditions

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  5                                     |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      T;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 0 0 1 0 0 0];

    internalField   uniform 278.15;     //278.15 K = 5 C

    boundaryField
    {
        left
        {
            type            zeroGradient;
        }
        right
        {
            type            zeroGradient;
        }
        top
        {
            //type            fixedValue;
            //value           uniform 273.15;
            type            inletOutlet;
            phi                     phi;
            inletValue      uniform 278.15;
        }
        bottom
        {
            type            codedFixedValue;
            value           uniform 873.15; 
            name            gaussShapeT;
            code            #{
                                scalarField x(this->patch().Cf().component(0)); 
                                double wGauss=200;
                                double x0=1000;
                                double Tmin=573;
                                double Tmax=873.15;
                                scalarField T(Tmin+(Tmax-Tmin)*exp(-(x-x0)*(x-x0)/(2*wGauss*wGauss)));
                                operator==(T);
                            #};
        }
        frontAndBack
        {
            type            empty;
        }
    }

    // ************************************************************************* //

The boundary conditions are again set for the patches that were defined in the blockMeshDict. Notice how the side are insulating (zeroGradient). The top has a  boundary condition called inletOutlet; it sets a constant inflow temperature (recharge of cold seawater) and assumes zeroGradient for the outflow (mimicing free fluid venting). The bottom boundary condition is special, it is set to codedFixedValue. The codedFixedValue BC allows "programming" a boundary condition on the fly. Here a gaussian-shapes constant temperature boundary condition is programmed. Note that :code:`x(this->patch().Cf().component(0))` is the x-coordinate of each FV face of the patch "bottom". 

Units are set by the dimensions keyword. The entries refer to the standard SI units [Kg m s K mol A cd]. By having a one in the fourth columns, the units of the defined properties has units of Kelvin.

We also need to set boundary conditions for pressure.

.. code-block:: bash 

    code 0/p

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 17, 33
    :name: lst:2dbox:bc_p
    :caption: Boundary conditions

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  5                                     |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      p;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [1 -1 -2 0 0 0 0];

    internalField   uniform 300e5;

    boundaryField
    {
        left
        {
            type            noFlux;
        }
        right
        {
            type            noFlux;
        }
        top
        {
            type            submarinePressure;
            value           uniform 300e5;
        }
        bottom
        {
            type            noFlux;
        }
        frontAndBack
        {
            type            empty;
        }
    }

    // ************************************************************************* //

The :code:`noFlux` boundary conditions, sets the pressure gradient to zero (horizontal direction) and hydrostatic (vertical direction), so that no flow occurs through these boundaries. The :code:`submarinePressure` boundary condition is provided by |foam| and sets the pressure according to water depth. Change it to fixedValue; we will discuss the special boundary conditions later.


Permeability
--------------------

In hydrothermal convection simulations, the fluid properties are given by the used EOS (details on this in the next lecture). What we need to set are the solid properties like permeability, solid density, solid specific heat, and porosity. These are set in two different files. Permeabilty is treated as a variable and is set in the 0 directory.

.. code-block:: bash 

    code 0/permeability

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 18
    :name: lst:2dbox:perm
    :caption: Permeability on hydrothermal flow simulations.

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  5.0                                   |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        location    "0";
        object      permeability;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 2 0 0 0 0 0];

    internalField   uniform 1e-14;

Again, check that you understand the units, which here add up to m^2. 


Case control
------------

Finally, we need to set some control parameters like the time step, run time, output writing. These kind of parameters are set in system/controlDict. Open it and explore the values. You will need to change the application to the new solver HydrothermalSinglePhaseDarcyFoam_xThermo.

.. code-block:: bash 

    code system/controlDict

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 16,17,20,21, 37
    :name: lst:2dbox:cdict
    :caption: controlDict of the Regular2DBox case.

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  5                                     |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      controlDict;
    }

    application HTFoam;
    startFrom startTime;
    startTime 0;
    stopAt endTime;
    endTime 86400000000;
    deltaT 864000;
    adjustTimeStep yes;
    maxCo           0.5;
    maxPorousCo     0.5; 
    maxDeltaT       86400000; 
    writeControl adjustableRunTime;
    writeInterval 864000000;
    purgeWrite 0;
    writeFormat ascii;
    writePrecision 6;
    writeCompression off;
    timeFormat general;
    timePrecision 14;
    runTimeModifiable true;
    libs 
    ( 
        "libHydrothermalBoundaryConditions.so"
    );



The solver we are using is called HydrothermalSinglePhaseDarcyFoam_xThermo. In addition, we are including the library "libHydrothermalBoundaryConditions.so", which provides special boundary conditions for submarine hydrothermal flow calculations.

.. warning::
    We have tuned the tutorials for runtime, so that you can run them on a standard laptop. However, the settings are not necessarily optimal for accuracy. For example, the time step is quite large and the maximum time step size if lmited by a Courant number computed for the Darcy velocity, while it should really be the pore velocity. If you want to do more accurate simulations, you should reduce the time step and the Courant number.



Solver controls
^^^^^^^^^^^^^^^
The numerical schemes and solver settings are set in the files fvSchemes and fvSolution, which are located in the system directory. Open them and check that you understand the settings. You can leave them as they are for now.

.. code-block:: foam 
    :linenos:
    :name: lst:2dbox:fvSolution
    :caption: fvSolution of the Regular2DBox case.


    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  10 (modern style)                      |
    |   \\  /    A nd           |                                                 |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      fvSolution;
    }

        // ------------------------------------------------------------------------- //
        // Linear solvers
        // ------------------------------------------------------------------------- //
        solvers
        {
            // Pressure
            p
            {
                solver          GAMG;
                tolerance       1e-10;       // linear solve absolute stop
                relTol          0.01;        // linear relative reduction (per solve)
                smoother        DICGaussSeidel;
            }

            pFinal
            {
                $p;
                tolerance       1e-10;       
                relTol          0;
            }

            // Temperature
            T
            {
                solver          GAMG;
                tolerance       1e-08;
                relTol          0.01;
                smoother        DILUGaussSeidel;
            }

            TFinal
            {
                $T;
                tolerance       1e-08;
                relTol          0;
            }
        }

        relaxationFactors
        {
            equations
            {
                p 1;
                T 1;
            }
        }

        PIMPLE
        {
            nOuterCorrectors         0;        
            nCorrectors              1;         
            nNonOrthogonalCorrectors 1;

            // Foundation OF expects single values here:
            residualControl
            {
                p 1e-4;
                T 1e-4;         
            }

            // (Optional) you can keep other PIMPLE keys here, but NOT dict-style RC.
        }

        PTCOUPLING
        {
            // Physical coupling criteria used by your init loop (custom)
            maxDeltaP         10;    // [Pa] absolute pressure change
            maxInitOuterIters 200;
            // number of iterations in main loop
            tightCouplingIters 3;    // K: number of mini Picard iterations per time step (0..3 typical)

        }




Running the case
----------------
Now we are finally ready to run our first test case. Just type this into your docker shell:

.. code-block:: bash 

    HTFoam

Notice how several directories are appearing, which contain the intermediate results. You can postprocess the case by simply opening the :code:`a.foam` file from paraview.


.. figure:: /_figures/RegularBox2D.*
   :align: center
   :name: fig:Regular2DBox_fig

   Results of the Regular2DBox example calculation.


.. warning::
    We have recently changed the solver settings and equation-of-state logic. Most of the tutorials in the cookbook directory are not yet updated and are missing, e.g., the new xThermoProperties dictionary. If you want to run these cases, you need to add this dictionary and set the parameters as shown above. We will update the cookbook in the near future.
