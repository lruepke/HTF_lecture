.. include:: /include.rst_

.. _L03_FirstCase:

Hydrothermal convection test case
==================================

Prepare case files
------------------
To get started we will run the **Regular2DBox** case from the cookbook directory of |foam|. This cookbooks describes how we can simulate a simple hydrothermal convection cell. It resolves hydrothermal convection driven by a gaussian-shaped constant temperature boundary condition at the bottom. 

Copy the  case into your shared working directory (probably $HOME/HydrothermalFoam_runs). You need to do this within the docker container (your right-hand shell in Visual Studio Code if you followed the recommended setup). Cd into your shared folder and type this:

.. code-block:: bash 
    :name: lst:cp2dBoxToWorkDir

    cd $HOME/HydrothermalFoam_runs
    cp -r $HOME/hydrothermalfoam-master/cookbooks/2d/Regular2DBox . 

Check out the directory structure shown in :numref:`lst:2dbox:tree`.

.. code-block:: bash 
    :linenos:
    :emphasize-lines: 3-5,10
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
 
The 0 directory now has entries for T (temperature) and p (pressure) our new primary variables, and for permeability, which we will discuss later. In addition, the constant directory has an entry thermophysicalProperties, which describes the solid properties.

.. tip::
    Most OpenFoam cases include scripts like :code:`run.sh` and :code:`clean.sh`. The :code:`run.sh` script is a good starting point for "understanding" a case. It lists all commands that have to be executed (e.g. meshing, setting of properties, etc.) to run a case. The :code:`clean.sh` script cleans up the case and deletes e.g. the mesh and all output directories. Have a look into these files and see if you understand them!

The 0 directory contains all initial and boundary conditions, the system folder contains all controlling parameter files, and the constant folder contains constant properties like the mesh - which we will create next. 

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


Transport properties
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


Next we look at the solid properties:

.. code-block:: bash 

    code constant/thermophysicalProperties

Check that you understand the units! Details can be found in the |foam| documentation.

Case control
------------

Finally, we need to set some control parameters like the time step, run time, output writing. These kind of parameters are set in system/controlDict. Open it and explore the values

.. code-block:: bash 

    code system/controlDict

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 16, 37-38
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

    application HydrothermalSinglePhaseDarcyFoam;
    startFrom latestTime;
    startTime 0;
    stopAt endTime;
    endTime 16912000000; //86400000000
    deltaT 864000;
    adjustTimeStep yes;
    maxCo           0.8; 
    maxDeltaT       86400000; 
    writeControl adjustableRunTime;
    writeInterval 86400000;
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
        "libHydroThermoPhysicalModels.so"
    );


The solver we are using is called HydrothermalSinglePhaseDarcyFoam. In addition, we are including two libraries "libHydrothermalBoundaryConditions.so"; these are part of |foam| and provide special boundary conditions for submarine hydrothermal flow calculations.

Running the case
----------------
Now we are finally ready to run our first test case. Just type this into your docker shell:

.. code-block:: bash 

    HydrothermalSinglePhaseDarcyFoam

Notice how several directories are appearing, which contain the intermediate results. You can postprocess the case by simply opening the :code:`a.foam` file from paraview.


.. figure:: /_figures/RegularBox2D.*
   :align: center
   :name: fig:Regular2DBox_fig

   Results of the Regular2DBox example calculation.

