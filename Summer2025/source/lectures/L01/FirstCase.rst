.. include:: /include.rst_

.. _L01_FirstCase:

A first test case - heat diffusion
==================================

.. admonition:: What do we need?

    - case files
    - a mesh
    - boundary/initial conditions
    - a solver
    - simulations controls (time step, run time, etc.)
    - post-processing

.. _theory_heat_diffusion:

Theory
------
The first case is the **flange** case from OpenFOAM tutorials, which illustrates how to solve for 3-D heat conduction in a complex-shaped solid. The equation we will be solving is the transient heat diffusion, a laplacian partial differential equation (PDE) that is here written for the simplfied case of constant properties (Eqn. :eq:`eq:laplacian_dif`). The equation is based on plugging Fourier's law (Eqn. :eq:`eq:fourier_law`) into a general energy conservation equation (Eqn. :eq:`eq:simple_energy_dif`)

.. math::
    :label: eq:fourier_law 
    
    \vec{q} = -\lambda_r \nabla T

.. math::
    :label: eq:simple_energy_dif 
    
    \rho C_{cp} \frac{\partial T}{\partial t}= -\nabla \cdot \vec{q}

.. math::
    :label: eq:laplacian_dif
    
    \frac{\partial T}{\partial t} = D \nabla^2 T   

Here  :math:`\rho` is density, :math:`C_{cp}` the specific heat,  :math:`\lambda_r` the heat conduction coefficient, and :math:`D=\frac{\lambda_r}{\rho C_{cp}}` the thermal diffusivity.

Prepare case files
------------------
Our first case is the **flange** case from the `OpenFoam tutorials <https://cfd.direct/openfoam/user-guide/>`_, which illustrates how to solve for 3-D heat conduction in a complex-shaped solid.

.. tip::
    It is good practice to never change the original tutorials, which reside in $FOAM_TUTORIALS but to always make a copy to your working directory.
 
To get started, just copy the whole case into your shared working directory (probably $HOME/HydrothermalFoam_runs). You need to do this within the docker container (your right-hand shell in Visual Studio Code if you followed the recommended setup). Cd into your shared folder and type this:

.. code-block:: bash 
    :name: lst:cpFlangeToWorkDir

    cd $HOME/HydrothermalFoam_runs
    cp -rf $FOAM_TUTORIALS/basic/laplacianFoam/flange . 

An OpenFoam case consist of a certain directory structure, which contains text files that control the case. The tree structure of flange case is shown in :numref:`lst:flange:tree`.

.. code-block:: bash 
    :linenos:
    :emphasize-lines: 8
    :name: lst:flange:tree
    :caption: File tree structure of the flange case.

    .
    ├── 0
    │   └── T
    ├── Allclean
    ├── Allrun
    ├── constant
    │   └── transportProperties
    ├── flange.ans
    └── system
        ├── controlDict
        ├── fvSchemes
        └── fvSolution

The 0 directory contains all initial and boundary conditions, the system folder contains all controlling parameter files, and the constant folder contains constant properties like the mesh - which we will create next. 

Mesh generation
---------------

The mesh of flange case is prepared as Ansys_ format (see line 8 of :numref:`lst:flange:tree`), you can use OpenFOAM mesh convert utility of :code:`ansysToFoam` to convert the ansys format to OpenFOAM format, the command is shown below,

.. code-block:: bash 

    ansysToFoam flange.ans -scale 0.001

The :code:`-scale` option is used to set scale factor of the whole mesh. 
Again, you can use command of :code:`tree` to explore what changes in the case folder, you can get something like :numref:`lst:flange_mesh:tree`, please notice the new files emphasized in the lines 7-13.

.. code-block:: bash 
    :linenos:
    :emphasize-lines: 7-13
    :name: lst:flange_mesh:tree
    :caption: File tree structure of the flange case.

    .
    ├── 0
    │   └── T
    ├── Allclean
    ├── Allrun
    ├── constant
    │   ├── polyMesh
    │   │   ├── boundary
    │   │   ├── faceZones
    │   │   ├── faces
    │   │   ├── neighbour
    │   │   ├── owner
    │   │   └── points
    │   └── transportProperties
    ├── flange.ans
    └── system
        ├── controlDict
        ├── fvSchemes
        └── fvSolution

After converting the mesh, you can use Paraview_ to visualize the mesh,

.. code-block:: bash

    touch a.foam
    paraview a.foam 

.. note:: 

    :code:`touch a.foam` means create a empty file named :code:`a.foam`. You might might have the :code:`touch` command in windows; then create the empty file inside the docker shell.

.. only:: html

    .. tab:: Mesh of flange

        .. figure:: /_figures/flange_mesh.png
            :align: center

            Mesh of the first case - flange.

    .. tab:: Interactive

        .. raw:: html
            
            <iframe src="../../_static/vtk_js/index.html?fileURL=flange_mesh.vtkjs" width="100%" height="500px"></iframe>
    
Boundary conditions
-------------------

Next we need to set boundary conditions. Open the file T inside the 0 directory from your local left-hand shell.

.. code-block:: bash 

    code 0/T

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 17,30,31
    :name: lst:flange_mesh:bc
    :caption: Boundary conditions

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
        object      T;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 0 0 1 0 0 0];

    internalField   uniform 273;

    boundaryField
    {
        patch1
        {
            type            zeroGradient;
        }

        patch2
        {
            type            fixedValue;
            value           uniform 273;
        }

        patch3
        {
            type            zeroGradient;
        }

        patch4
        {
            type            fixedValue;
            value           uniform 573;
        }
    }

    // ************************************************************************* //

Notice how boundary conditions are set to parts of the mesh that have an identifier like patch2. Those identifiers are defined during the meshing procedure. Here they were defined in the ansys file. 

Units are set by the dimensions keyword. The entries refer to the standard SI units [Kg m s K mol A cd]. By having a "1" in the fourth columns, the units of the defined properties has units of Kelvin.

Transport properties
--------------------

Next we need to set the value of the thermal diffusivity. 

.. code-block:: bash 

    code constant/transportProperties

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 18
    :name: lst:flange_mesh:tp
    :caption: Transport properties

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
        class       dictionary;
        location    "constant";
        object      transportProperties;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    DT              DT [0 2 -1 0 0 0 0] 4e-05;


    // ************************************************************************* //

Again, check that you understand the units, which here add up to m^2/s. 

Case control
------------

Finally, we need to set some control parameters like the time step, run time, output writing. These kind of parameters are set in system/controlDict. Open it and explore the values

.. code-block:: bash 

    code system/controlDict

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 18
    :name: lst:flange_mesh:cd
    :caption: controlDict of the flange case.

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
        class       dictionary;
        location    "system";
        object      controlDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    application     laplacianFoam;

    startFrom       latestTime;

    startTime       0;

    stopAt          endTime;

    endTime         3;

    deltaT          0.005;

    writeControl    runTime;

    writeInterval   0.1;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;

You can find detailed descriptions of the different parameters in the `OpenFoam User Guide <https://cfd.direct/openfoam/user-guide/>`_. Note that the solver we are using is called laplacianFoam. 

Running the case
----------------
Now we are finally ready to run our first test case. Just type this into your docker shell:

.. code-block:: bash 

    laplacianFoam

Notice how several directories are appearing, which contain the intermediate results. You can postprocess the case by simply opening the :code:`a.foam` file from Paraview.

