Navier-Stokes flow with OpenFoam
================================

In this lecture we will solve the famous Navier-Stokes equation for flow. At this point, we will not worry about the detailed equations we are solving. Let's just accept that the Navier-Stokes equation(s) allow us to solve for the motion of a fluid subject to boundary conditions (simply speaking). 


2D cavity flow
------------------------
Our first example of simulating Navier-Stokes flow is flow within a 2D cavity. Flow is driven by a kinematic boundary condition at one side, while all other sides are walls with zero velocity. The model setup is shown in :numref:`fig:2dcavity_setup` .

.. figure:: /_figures/2d_cavity_setup.*
   :align: center
   :name: fig:2dcavity_setup

   Two-dimensional cavity model geometry and boundary conditions.


Preparing the case
^^^^^^^^^^^^^^^^^^

Two-dimensional cavity flow is included in the official tutorials of OpenFoam. We get started by copying the tutorial case into our work directory. You need to do this from your docker shell.

.. code-block:: bash 
    :name: lst:cp2dCavityToWorkDir

    cd $HOME/HydrothermalFoam_runs
    cp $FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity ./cavity2D

Check the directory structure:

.. code-block:: bash 
    :name: lst:2dcavity_dir
    :emphasize-lines: 2,7,9

    .
    ├── 0
    │   ├── U
    │   └── p
    ├── a.foam
    ├── clean.sh
    ├── constant
    │   └── transportProperties
    └── system
        ├── blockMeshDict
        ├── controlDict
        ├── fvSchemes
        └── fvSolution

The 0 directory now contains the boundary and initial conditions for pressure and velocity (our primary variables), the constant directory has the transport properties (here viscosity), and the system directory contains information on the mesh, time steps, numerical schemes etc.

The two files :code:`run.sh` and :code:`clean.sh` are actually not included and we need to create them. These files are usually part of every OpenFoam case; :code:`run.sh` includes all commands that need to executed to run the case and :code:`clean.sh` cleans the case (removes mesh and output directories).

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

Make the scripts executable.

.. code-block:: bash 
    :name: lst:2dCavitychmod

    chmod u+x clean.sh run.sh

Making the mesh
^^^^^^^^^^^^^^^
We will use OpenFoam's blockMesh utility to make a simple 2D mesh. The corresponding :code:`blockMeshDict` file that has all the meshing information is located in the system folder.

.. figure:: /_figures/cavity2d_bm.*
   :align: center
   :name: fig:cavity2d_bm

   Structure of the blockMeshDict

First we need to define the vertices of the mesh, the nodes.

.. figure:: /_figures/cavity2d_vertices.*
   :align: center
   :name: fig:cavity2d_vertices

   Numbering of the vertices.

The next step is define the connectivity between the vertices in order to describe the modeling domain.

.. figure:: /_figures/cavity2d_vorder.*
   :align: center
   :name: fig:cavity2d_vorder

   The order by which the vertices are passed to the hex command matters!

.. admonition:: Order of vertices

    The OpenFoam documentation provides a nice description of the vertices ordering.
    
    * the axis origin is the first entry in the block definition, vertex 0 in our example
    * the x direction is described by moving from vertex 0 to vertex 1
    * the y direction is described by moving from vertex 1 to vertex 2
    * vertices 0, 1, 2, 3 define the plane z = 0
    * vertex 4 is found by moving from vertex 0 in the z direction
    * vertices 5,6 and 7 are similarly found by moving in the z direction from vertices 1,2 and 3 respectively.

Next boundary patches are defined and labeled in the blockMeshDict. Also here care must be take to provide the vertices in a consistent order (right-hand coordinate system). Two easy ways to remember this is to:

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

You can visualize the mesh using paraview

.. code-block:: bash 
    :name: lst:2dCavity_vizm

    touch a.foam 
    paraview a.foam 

Boundary conditions
^^^^^^^^^^^^^^^^^^^

We now have velocity and pressure as primary variables and need to set initial and boundary conditions for them. First we look at the velocity boundary conditions:

.. code-block:: bash 

    code 0/u 


.. figure:: /_figures/cavity2d_u.*
   :align: center
   :name: fig:cavity2d_u_fig

   Velocity boundary conditions. The front and back sides are set to empty because we are doing a 2D calculation.

Next we look into the pressure boundary conditions.

.. code-block:: bash 

    code 0/p 

.. code-block:: foam 
    :name: lst:2dcavity_p
    :emphasize-lines: 17
    :linenos:
    :caption: Pressure boundary conditions. Front and back are of type "emtpy" for 2-D runs. 

    /*--------------------------------*- C++ -*----------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     | Website:  https://openfoam.org
      \\  /    A nd           | Version:  7
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
    One has to be careful about the dimensions of pressure in OpenFoam. In incompressible runs, like we are doing here, the pressure is usually the relative pressure :math:`\frac{p}{\rho}` and has units :math:`\frac{m^2}{s^2}` 


Run controls
^^^^^^^^^^^^^^^^^^^

The time stepping, run time, and output frequency are again set in :code:`system/controlDict`. Open it and check that you understand the entires. 

In case you wondered how OpenFoam is solving the equations. We will cover the details later in the course, but you can have a preview by opening the :code:`system/fvSchemes` file. In this dictionary, the various discretization schemes can be set. :numref:`fig:cavity2d_num_fig` gives some further explanations.


.. figure:: /_figures/cavity2_num.*
   :align: center
   :name: fig:cavity2d_num_fig

   The exact discretization schemes can be set in :code:`system/fvSchemes`.


Time to run the case! Just start the solver

.. code-block:: bash 

    icoFoam

.. admonition:: Technical detail

   The :code:`icoFoam` solver is actually using the PISO algorithm and the not the SIMPLE algorithm described above. SIMPLE is a steady-state solver, while PISO can also resolve transient flows. PISO also make some smart corrections for better convergence. We will learn about these technical issues later in this course. 


Visualization
^^^^^^^^^^^^^^^^^^^
Open paraview and look at the results.