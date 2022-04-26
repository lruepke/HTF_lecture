.. include:: /include.rst_

.. _L03_FirstCase:

Computing effective permeability
==================================

Objective
---------

Starting from a digital representation of a sample's pore space (typically an image produced by a cT-scan), we want to compute its effective permeability. Or put differently, we will make a **direct** simulation of flow on the pore level and post-process it for extracting the effective permeability for simplified **continuum** simulations using Darcy's law.


.. figure:: /_figures/porousModel.*
   :align: center
   :name: fig:porousModel_fig
   :figwidth: 70%

   Synthetic image of the pore geometry. Pores are white; solid is black.


To compute the effective permeability, we will apply constant pressure boundary conditions and evaluate the flow rate through the sample. Once we have that, we can re-arrange Darcy's law to solve for permeability:

.. math::
    :label: eq:darcy_perm 

        K_{i,j} \, = \, \mu \Biggl( \frac{ \Delta P}{ \Delta x_j} \Biggr)^{-1} \, \Biggl( \frac{1}{V} \int_V u_i \, dV \Biggr) \; .


Workflow
--------

Ok, let's do it!

Mesh generation
^^^^^^^^^^^^^^^
The first major step is to generate a mesh of the pores space shown in :numref:`fig:porousModel_fig`. For this purpose, we will use `OpenFOAM's snappyHexMesh tool <https://cfd.direct/openfoam/user-guide/v7-snappyHexMesh/#x27-1970005.4>`_ . It allows meshing arbitrary geometries and is very powerful. Unfortunately, it can also be infuriating to use as it asks for many user-defined parameters and can be quite picky about the choices made. 

.. tip::
    We will not go into the details of snappyHexMesh (SHM). If you want to use it and/or understand it, a good starting point is the user-guide linked above. Another great resource are the `Rock Vapor Classic tutorial series <https://holzmann-cfd.com/community/training-videos/openfoam-usage/rock-vapor-classic>`_ .

In a nutshell, snappyHexMesh (SHM) is about starting from a blockMesh (as in the previous lecture) and then cutting-out the solid grain (described by a triangulated surface), and then *snapping* the mesh to this surface. A typical way to describe this surface is an .stl file - a typical file format for triangulated surfaces that is also often used for 3D printing.


.. figure:: /_figures/figure_workflow.*
   :align: center
   :name: fig:figure_workflow_fig
   :figwidth: 85%

   Workflow illustrating the meshing process.


The steps involved are shown in :numref:`fig:figure_workflow_fig` . Starting from an image, an stl file is created that is then used during the meshing process. Most of the steps will rely on paraview filters.

    #. start with an image (A).
    #. save it as a .vti file that is easily understood by Paraview. We use `porespy <https://porespy.org/>`_ for this step.
    #. load the vti file into paraview and use the *clip* (B) and *triangulation* (B) filers to created a surface of the pore space (C).
    #. save this surface as a stl file


Python pre-processing
^^^^^^^^^^^^^^^^^^^^^

Let's work through the steps involved and assume we received a 2-D image of scanned pore space ( :numref:`fig:figure_workflow_fig` A). We need to translate it into something that Paraview understands, so that we can do the segmentation and surface generation. We will use porespy for it and the first steps are to install porespy into our python virtual environment (we should already have PIL, which is also needed).

.. code-block:: bash

      conda activate "your_environment_name"
      conda install -c conda-forge porespy

.. tip::

    If conda install fails, you can also use :code:`pip install porespy`
    
    
Now we are good to go and it's time to download the data. The complete openFOAM case can be downloaded from :download:`here <cases/DRP_permeability_2D.zip>` . Next we import the .png file from the :code:`geometry` folder and convert it to a .vti file, which is the `Visualization Toolkit's <https://vtk.org>`_ format for storing image data. Note that vti files can also store series of images, which is important when doing this in 3-D.

.. code-block:: python
    
    import numpy as np   
    from PIL import Image
    import porespy as ps

    impath = 'geometry/'         # image path
    imname = 'porousModel.png'   # file name 

    image = Image.open('%s/%s'%(impath,imname)).convert("L")                   # open image
    arr = np.asarray(image)                                                    # convert image to array
    ps.io.to_vtk(np.array(arr, dtype=int)[:, :, np.newaxis], 'porous_model')   # use porespy and save to .vti format

You can put this little script into a jupyter notebook, or save it as .py file to be run from the command line. After running it, you should  have a file :code:`porous_model.vti` in the folder where you executed the script. Now comes the segementation and triangulation part. We can use paraview's python interface for this (or do everything by hand using the graphical user interface).

.. code-block:: python

    # workflow as python code using the paraview.simple module
    from paraview.simple import *
    
    def write_stl(vti_file, stl_file):
        data = OpenDataFile('%s.vti'%vti_file)
        clip1 = Clip(data, ClipType = 'Scalar', Scalars = ['CELLS', 'im'], Value = 127.5, Invert = 1) 
        extractSurface1 = ExtractSurface(clip1)
        triangulate1 = Triangulate(extractSurface1)
        SaveData(stl_file, proxy = triangulate1)
    
    vti_file = 'porous_model'      # input .vti file
    stl_file = 'porous_model.stl'  # output .stl file

    write_stl(vti_file, stl_file)


Getting the paraview.simple module to work in a virtual environment is a challenge. The easiest workaround is to use paraview's python shell and execute the script there:

.. figure:: /_figures/paraview_python.*
   :align: center
   :name: fig:figure_workflow_fig
   :figwidth: 85%

   Using the paraview python shell.


An easy way to do this is to use these commands:

.. code:: python

    import os 

    os.chdir("your_case_directory")
    exec(open("your_file_name.py").read())
    

This will create a .stl file, which we will use in the meshing process. Let's do some clean up and move the vti file into :code:`./geometry` and the stl file into :code:`./constant/triSurface`, where openFOAM expects it. This assumes that you are in the case directory.

.. code:: bash

    mv ./porous_model.vti ./geometry/
    mv ./porous_model.stl ./constant/triSurface/


OpenFOAM case
^^^^^^^^^^^^^

Great, back to openfoam for the final mesh making! Making the mesh with SHM is a two-step process. First we make a standard blockMesh background mesh. This is, as usual, controlled by :code:`system/blockMeshDict`:

.. code-block:: foam 
    :caption: blockMeshDict
    :emphasize-lines: 24, 25, 40, 49
    :linenos:


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
        object      blockMeshDict;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    convertToMeters 1;

    lx0 0;
    ly0 0;
    lz0 0;

    lx1 1196;
    ly1 1494;
    lz1 1;

    vertices
    (
        ($lx0 $ly0 $lz0)   //0
        ($lx1 $ly0 $lz0)   //1
        ($lx1 $ly1 $lz0)   //2
        ($lx0 $ly1 $lz0)   //3
        ($lx0 $ly0 $lz1)   //4
        ($lx1 $ly0 $lz1)   //5
        ($lx1 $ly1 $lz1)   //6
        ($lx0 $ly1 $lz1)   //7
    );

    blocks
    (
        hex (0 1 2 3 4 5 6 7) (315 390 1) simpleGrading (1 1 1)
    );

    edges
    (
    );

    boundary
    (
        top
        {
            type symmetryPlane;
            faces
            (
                (7 6 3 2)
            );
        }

        inlet
        {
            type wall;
            faces
            (
                (0 4 7 3)
            );
        }

        bottom
        {
            type symmetryPlane;
            faces
            (
                (1 5 4 0)
            );
        }

        outlet
        {
            type patch;
            faces
            (
                (1 2 6 5)
            );
        }


        frontAndBack
        {
            type empty;
            faces
            (
                (0 3 2 1)
                (4 5 6 7)
            );
        }
    );


Notice how we set the vertical and horizontal extents to 1196 and 1494, which is the pixel resolution of the image. We will scale it later to physical dimensions. Notice also how the box is described and how boundary conditions are applied.

Now have a look at :code:`system/snappyHexMeshDict`. We will not go into details here, just explore the general structure yourself if you are interested using the resources linked above.

Now it's time to make the mesh. Run each of the steps below individually and check out the results in paraview!

.. code:: bash

    blockMesh
    snappyHexMesh -overwrite
    checkMesh -allTopology -allGeometry
    transformPoints -scale "(1e-6 1e-6 1e-6)"

The final conversion turns everything in micrometer (:math:`10^{-6} m`). Check out the final mesh in paraview, it should look like this:




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

