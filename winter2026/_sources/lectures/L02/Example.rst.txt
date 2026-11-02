Navier-Stokes flow with OpenFOAM
================================

In this lecture we will solve the famous Navier-Stokes equation for flow. There are many different forms of the Navier-Stokes equations, depending on the assumptions made. Here we solve the form that is treated in the ``icoFoam`` solver in OpenFOAM. It treats **incompressible**, **laminar**, and **transient** fluid flow and solves the following equations:

.. math::

    \frac{\partial \mathbf{U}}{\partial t} + (\mathbf{U} \cdot \nabla)\mathbf{U} = -\nabla p + \nu \nabla^2 \mathbf{U}

.. math::

    \nabla \cdot \mathbf{U} = 0

where:

- :math:`\mathbf{U}` is the velocity vector field,
- :math:`p` is the kinematic pressure (i.e., pressure divided by density),
- :math:`\nu` is the kinematic viscosity,
- :math:`t` is time.

The first equation is the **momentum conservation equation**:

- :math:`\partial \mathbf{U} / \partial t` represents unsteady acceleration (how velocity changes over time),
- :math:`(\mathbf{U} \cdot \nabla)\mathbf{U}` is the advection (fluid moving itself),
- :math:`-\nabla p` is the pressure gradient force,
- :math:`\nu \nabla^2 \mathbf{U}` is the viscous diffusion of momentum.

The second equation enforces **incompressibility**, meaning the volume of fluid elements does not change over time.

In ``icoFoam``, these equations are solved using the **finite volume method** and the **PISO algorithm** (Pressure Implicit with Splitting of Operators) to ensure pressure–velocity coupling and divergence-free velocity fields at each time step. We will get into the details of the numerical methods later in the course. For now, we will just use the solver and its capabilities.


2D cavity flow
------------------------
Our first example of simulating Navier-Stokes flow is flow within a 2D cavity. Flow is driven by a kinematic boundary condition at one side, while all other sides are walls with zero velocity. The model setup is shown in :numref:`fig:2dcavity_setup` .

.. figure:: /_figures/2d_cavity_setup.*
   :align: center
   :name: fig:2dcavity_setup

   Two-dimensional cavity model geometry and boundary conditions. Note a mistake in the figure: the domain size is actually 0.1 m by 0.1 m in the tutorial.


Preparing the case
^^^^^^^^^^^^^^^^^^

Two-dimensional cavity flow is included in the official tutorials of OpenFOAM. We get started by copying the tutorial case into our work directory. You need to do this in a terminal in the container (e.g. in Visual Studio Code).

.. code-block:: bash 
    :name: lst:cp2dCavityToWorkDir

    cd $HOME/HydrothermalFoam_runs
    cp -r $FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity ./cavity2D

Go into the case and check the directory structure with :code:`tree`:

.. code-block:: bash

    cd cavity2D
    tree

.. code-block:: bash
    :name: lst:2dcavity_dir
    :emphasize-lines: 2,5,7

    .
    ├── 0
    │   ├── U
    │   └── p
    ├── constant
    │   └── physicalProperties
    └── system
        ├── blockMeshDict
        ├── controlDict
        ├── fvSchemes
        └── fvSolution

The 0 directory now contains the boundary and initial conditions for pressure and velocity (our primary variables), the constant directory has the physical properties (here the kinematic viscosity), and the system directory contains information on the mesh, time steps, numerical schemes etc.

The two files :code:`run.sh` and :code:`clean.sh` are actually not included and we need to create them. These files are usually part of every OpenFOAM case; :code:`run.sh` includes all commands that need to be executed to run the case and :code:`clean.sh` cleans the case (removes mesh and output directories).

To create them, make sure you are in the case directory and open a new, empty file in the editor of Visual Studio Code:

.. code-block:: bash

    cd $HOME/HydrothermalFoam_runs/cavity2D
    code run.sh

Copy the content of the first listing below (:code:`run.sh`) into the editor and save the file (CTRL/CMD + S). Then do the same for :code:`clean.sh` with the second listing:

.. code-block:: bash

    code clean.sh

.. code-block:: bash 
    :linenos:
    :caption: The run.sh file.

    #!/bin/sh
    cd ${0%/*} || exit 1    # Run from this directory

    # Source tutorial run functions
    . $WM_PROJECT_DIR/bin/tools/RunFunctions

    application=`getApplication`

    ./clean.sh
    runApplication blockMesh
    runApplication $application


.. code-block:: bash 
    :linenos:
    :name: lst:2dcav_clean:tree
    :caption: The clean.sh file.

    #!/bin/sh
    cd ${0%/*} || exit 1 # run from this directory

    # Source tutorial run functions
    . $WM_PROJECT_DIR/bin/tools/CleanFunctions

    # Clean the case
    cleanCase

Make the scripts executable.

.. code-block:: bash 
    :name: lst:2dCavitychmod

    chmod u+x clean.sh run.sh

Making the mesh
^^^^^^^^^^^^^^^
We will use OpenFOAM's blockMesh utility to make a simple 2D mesh. The corresponding :code:`blockMeshDict` file that has all the meshing information is located in the system folder.

.. figure:: /_figures/cavity2d_bm.*
   :align: center
   :name: fig:cavity2d_bm

   Structure of the blockMeshDict

First we need to define the vertices of the mesh, the nodes.

.. figure:: /_figures/cavity2d_vertices.*
   :align: center
   :name: fig:cavity2d_vertices

   Numbering of the vertices.

The next step is to define the connectivity between the vertices in order to describe the modeling domain.

.. figure:: /_figures/cavity2d_vorder.*
   :align: center
   :name: fig:cavity2d_vorder

   The order by which the vertices are passed to the hex command matters!

.. admonition:: Order of vertices

    The OpenFOAM documentation provides a nice description of the vertices ordering.
    
    * the axis origin is the first entry in the block definition, vertex 0 in our example
    * the x direction is described by moving from vertex 0 to vertex 1
    * the y direction is described by moving from vertex 1 to vertex 2
    * vertices 0, 1, 2, 3 define the plane z = 0
    * vertex 4 is found by moving from vertex 0 in the z direction
    * vertices 5,6 and 7 are similarly found by moving in the z direction from vertices 1,2 and 3 respectively.

Next boundary patches are defined and labeled in the blockMeshDict. Also here care must be taken to provide the vertices in a consistent order (right-hand coordinate system). Two easy ways to remember this are:

    * apply the right-hand rule, which means if the thumb of your right hand points to the outside of a face, the numbering has to follow your fingers.

    * or, looking onto a face and starting from any vertex, the numbering has to be counter-clockwise.

.. figure:: /_figures/cavity2d_bounds.*
   :align: center
   :name: fig:cavity2d_bounds

   Assigning boundary labels and types.

Now we are ready to run the :code:`blockMesh` utility and create the mesh

.. code-block:: bash 
    :name: lst:2dCavityrbm

    blockMesh

You can visualize the mesh in ParaView. As in the first lecture, create the dummy file :code:`a.foam` in the case directory:

.. code-block:: bash
    :name: lst:2dCavity_vizm

    touch a.foam

Then open ParaView on your computer and load :code:`a.foam` from your course folder under :code:`HydrothermalFoam_runs/cavity2D`.

Boundary conditions
^^^^^^^^^^^^^^^^^^^

We now have velocity and pressure as primary variables and need to set initial and boundary conditions for them. First we look at the velocity boundary conditions:

.. code-block:: bash 

    code 0/U

.. code-block:: foam
    :name: lst:2dcavity_U
    :emphasize-lines: 16,18,24-25,30,35
    :linenos:
    :caption: Velocity boundary conditions. The front and back sides are set to "empty" because we are doing a 2-D calculation.

    /*--------------------------------*- C++ -*----------------------------------*\
      =========                 |
      \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
       \\    /   O peration     | Website:  https://openfoam.org
        \\  /    A nd           | Version:  10
         \\/     M anipulation  |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        format      ascii;
        class       volVectorField;
        object      U;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 1 -1 0 0 0 0];

    internalField   uniform (0 0 0);

    boundaryField
    {
        movingWall
        {
            type            fixedValue;
            value           uniform (1 0 0);
        }

        fixedWalls
        {
            type            noSlip;
        }

        frontAndBack
        {
            type            empty;
        }
    }

    // ************************************************************************* //

The dimensions of the field are m/s (line 16). Velocity is a vector with three components, so the initial condition :code:`uniform (0 0 0)` (line 18) means that the fluid is at rest everywhere in the domain. The lid (:code:`movingWall`) moves with 1 m/s in the x direction, and :code:`noSlip` sets the velocity to zero at the fixed walls.

Next we look into the pressure boundary conditions.

.. code-block:: bash 

    code 0/p 

.. code-block:: foam 
    :name: lst:2dcavity_p
    :emphasize-lines: 17
    :linenos:
    :caption: Pressure boundary conditions. Front and back are of type "empty" for 2-D runs. 

    /*--------------------------------*- C++ -*----------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     | Website:  https://openfoam.org
      \\  /    A nd           | Version:  10
       \\/     M anipulation  |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      p;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 2 -2 0 0 0 0];

    internalField   uniform 0;

    boundaryField
    {
        movingWall
        {
            type            zeroGradient;
        }

        fixedWalls
        {
            type            zeroGradient;
        }

        frontAndBack
        {
            type            empty;
        }
    }

    // ************************************************************************* //

.. tip::
    One has to be careful about the dimensions of pressure in OpenFOAM. In incompressible runs, like we are doing here, the solver uses the kinematic pressure :math:`p/\rho` (the pressure divided by the density), which has units of m²/s². This is why the dimensions in :code:`0/p` are :code:`[0 2 -2 0 0 0 0]` and not those of Pa.


Run controls
^^^^^^^^^^^^^^^^^^^

The time stepping, run time, and output frequency are again set in :code:`system/controlDict`. Open it and check that you understand the entries. 

.. code-block:: foam 
    :name: lst:2dcavity_controlDict
    :linenos:
    :caption: controlDict of the cavity flow tutorial. 

    /*--------------------------------*- C++ -*----------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     | Website:  https://openfoam.org
      \\  /    A nd           | Version:  10
       \\/     M anipulation  |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        format      ascii;
        class       dictionary;
        location    "system";
        object      controlDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    application     icoFoam;

    startFrom       startTime;

    startTime       0;

    stopAt          endTime;

    endTime         0.5;

    deltaT          0.005;

    writeControl    timeStep;

    writeInterval   20;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;


    // ************************************************************************* //

Here is a short explanation of the most important entries:


- ``application``: Solver used (here, ``icoFoam`` for incompressible laminar flow).
- ``startFrom``: Start method (``startTime`` = use value from ``startTime``).
- ``startTime``: Initial simulation time.
- ``stopAt``: Stop method (``endTime`` = use value from ``endTime``).
- ``endTime``: Final simulation time.
- ``deltaT``: Time step size.
- ``writeControl``: When to write output (``timeStep`` = every N time steps).
- ``writeInterval``: Write output every N time steps.
- ``purgeWrite``: Max number of saved time steps (``0`` = keep all).
- ``writeFormat``: Output file format (``ascii`` or ``binary``).
- ``writePrecision``: Number of digits in numerical output.
- ``writeCompression``: Compress output files (``off`` or ``on``).
- ``timeFormat``: Time label format (e.g., ``general``, ``fixed``).
- ``timePrecision``: Precision of time labels.
- ``runTimeModifiable``: Allow runtime dictionary edits (``true`` = yes).


These are the most important entries. There are many more options available, which can be found in the `OpenFOAM User Guide <https://doc.cfd.direct/openfoam/user-guide-v10/controldict>`_


In case you wondered how OpenFOAM is solving the equations. We will cover the details later in the course, but you can have a preview by opening the :code:`system/fvSchemes` file. In this dictionary, the various discretization schemes can be set. :numref:`fig:cavity2d_num_fig` gives some further explanations.


.. figure:: /_figures/cavity2_num.*
   :align: center
   :name: fig:cavity2d_num_fig

   The exact discretization schemes can be set in :code:`system/fvSchemes`.


Time to run the case! Use the script you have written, which creates the mesh and runs the solver:

.. code-block:: bash

    ./run.sh

The output of the solver is written to the file :code:`log.icoFoam`. Have a look at it, e.g. with :code:`code log.icoFoam`.

Visualization
^^^^^^^^^^^^^^^^^^^
Open the case in ParaView and look at the results. They are actually quite boring because the flow is steady and laminar; and because we have not first thought about the expected results. Let's do that better in the next exercise.

