.. include:: /include.rst_

.. _L03_FirstCase:

Flow on the pore scale
===========================

Objective
---------

Starting from a digital representation of a sample's pore space (typically an image produced by a cT-scan), we want to compute the permeability of the sample. Put differently, we will make a **direct** simulation of flow on the pore level and post-process it for extracting the permeability to be used, for example, in simplified **continuum** simulations using Darcy's law.

.. figure:: /_figures/porousModel.*
   :align: center
   :name: fig:porousModel_fig
   :figwidth: 70%

   Synthetic image of the pore geometry. Pores are white; solids are black.


To compute permeability, we will apply constant pressure boundary conditions and evaluate the flow rate through the sample. Once we have that, we can re-arrange Darcy's law to solve for permeability:

.. math::
    :label: eq:darcy_perm 

        K_{i,j} \, = \, \mu \Biggl( \frac{ \Delta P}{ \Delta x_j} \Biggr)^{-1} \, \Biggl( \frac{1}{V} \int_V u_i \, dV \Biggr) \; .


Ok, let's do it!

Mesh generation
---------------

The first major step is to generate a mesh of the pore space (:numref:`fig:porousModel_fig`). For this purpose, we will use `OpenFOAM's snappyHexMesh tool <https://cfd.direct/openfoam/user-guide/v7-snappyHexMesh/#x27-1970005.4>`_. It allows meshing arbitrary and complex geometries and is very powerful. Unfortunately, it can also be infuriating to use as it asks for many user-defined parameters and can be quite picky about the choices made. 

.. tip::
    We will not go into the details of snappyHexMesh (SHM) works. If you want to use and/or understand it, a good starting point is the user guide linked above. Another great resource is the `Rock Vapor Classic tutorial series <https://holzmann-cfd.com/community/training-videos/openfoam-usage/rock-vapor-classic>`_.

In a nutshell, snappyHexMesh (SHM) is about starting from a blockMesh (as in the previous lecture), cutting out the solid grain (described by a triangulated surface), and then *snapping* the mesh to this surface. A typical way to describe this surface is an .stl file - a file format for triangulated surfaces that is often used for 3D printing.


.. figure:: /_figures/figure_workflow.*
   :align: center
   :name: fig:figure_workflow_fig
   :figwidth: 85%

   Workflow illustrating the meshing process.


The steps involved are shown in :numref:`fig:figure_workflow_fig` . Starting from an image, an stl file is created that is then used during the meshing process. Most of the steps will rely on paraview filters and the workflow is this.

    #. Start with an image (A).
    #. Save it as a .vti file that is easily understood by Paraview. We use `porespy <https://porespy.org/>`_ for this step.
    #. Load the vti file into paraview and use the *clip* (B) and *triangulation* (B) filers to created a surface of the pore space (C).
    #. Save this surface as a stl file


Python pre-processing
^^^^^^^^^^^^^^^^^^^^^

Let's work through the steps involved and assume we received a 2-D image of scanned pore space ( :numref:`fig:figure_workflow_fig` A). We need to translate it into something that Paraview understands, so that we can do the segmentation and surface generation. We will use porespy for it and the first step is to install porespy into our python virtual environment (we should already have PIL, which is also needed).

.. code-block:: bash

      conda activate "your_environment_name"
      conda install -c conda-forge porespy

.. tip::

    If conda install fails, you can also use :code:`pip install porespy`
    
    
Now we are good to go and it's time to download the data. The complete openFOAM case can be downloaded from :download:`here <cases/DRP_permeability_2D.zip>` . 

Copy the case into your shared working directory (probably $HOME/HydrothermalFoam_runs). You need to do this within the docker container (your right-hand shell in Visual Studio Code if you followed the recommended setup).

Check out the directory structure shown in :numref:`lst:drp_case:tree`.

.. code-block:: bash 
    :name: lst:drp_case:tree
    :caption: File tree structure of the Digital Rock Physics case.

    .
    ├── 0
    │   ├── U
    │   └── p
    ├── a.foam
    ├── clean.sh
    ├── constant
    │   ├── transportProperties
    │   ├── triSurface
    │   └── turbulenceProperties
    ├── geometry
    │   └── porousModel.png
    ├── run.sh
    └── system
        ├── blockMeshDict
        ├── controlDict
        ├── fvSchemes
        ├── fvSolution
        ├── meshQualityDict
        └── snappyHexMeshDict



.. tip::
    Most OpenFoam cases include scripts like :code:`run.sh` and :code:`clean.sh`. The :code:`run.sh` script is a good starting point for "understanding" a case. It lists all commands that have to be executed (e.g. meshing, setting of properties, etc.) to run a case. The :code:`clean.sh` script cleans up the case and deletes e.g. the mesh and all output directories. Have a look into these files and see if you understand them!


Next we import the .png file from the :code:`geometry` folder and convert it to a .vti file, which is the `Visualization Toolkit's <https://vtk.org>`_ format for storing image data. Note that vti files can also store series of images, which is important when doing this in 3-D.

.. code-block:: python
    
    import numpy as np
    from PIL import Image
    import porespy as ps
    import os

    impath = 'geometry/'         # image path
    imname = 'porousModel.png'   # file name

    # 1. use PIL to load the image
    # convert("L") converts RGB into greyscale, L = R * 299/1000 + G * 587/1000 + B * 114/1000
    image = Image.open(os.path.join(impath, imname)).convert("L")

    # 2. convert image into numpy array with "normal" integer values
    arr = np.asarray(image, dtype=int)

    # 3. save as vti (VTK's image format)
    # 3.1 make 3D by repeating in new axis=2 dimension (openfoam and porespy want 3D)
    arr_stacked = np.stack((arr,arr), axis=2)

    # 3.2 save as vti using porespy                                              
    ps.io.to_vtk(arr_stacked, 'geometry/porous_model') 

You can put this little script into a jupyter notebook or save it as .py (e.g. png2vti.py), then activate the right kernel (e.g. :code:`conda activate py3_htf_class`), and do this:

.. code-block:: bash

    python png2vti.py


After running it, you should  have a file :code:`porous_model.vti` in the geometry folder. 

.. tip::
    If you are interested in the vti file format, this is what we have just written to porous_model.vti

    .. code-block:: bash

        <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt64">
            <ImageData WholeExtent="0 1196 0 1494 0 2" Origin="0 0 0" Spacing="1 1 1" Direction="1 0 0 0 1 0 0 0 1">
                <Piece Extent="0 1196 0 1494 0 2">
                <PointData>
                </PointData>
                <CellData Scalars="im">
                <DataArray type="Int64" Name="im" format="ascii" RangeMin="0" RangeMax="255">
                    255 255 255 255 255 255
                    255 255 255 255 255 255
                    255 255 255 255 255 255
                    ...

    Notice how this is written as cell data. The total extents are 1197x1495x3, the voxel size is 1, and the cell data is 1196x1494x2 (our two images).


Now comes the segmentation and triangulation part to make an stl file that openFOAM understands. We will use paraview for this and there are different ways of doing this:

    #. use paraview's graphical user interface
    #. use paraview's inbuild python shell
    #. install paraview into a conda environment and use a jupyter notebook

.. figure:: /_figures/paraview_python.*
   :align: center
   :name: fig:paraview_figure_fig
   :figwidth: 85%

   Using the paraview python shell.


Here is some python code that uses the paraview python bindings. You can save it into your case directory. Or, you use the paraview GUI to create the stl file.

.. code-block:: python

    # workflow as python code using the paraview.simple module
    from paraview.simple import *
    import os

    def write_stl(vti_file, stl_file):
        # 1. load vti file
        data            = OpenDataFile(vti_file)
        # 2. clip at some intermediate value (we have 0 and 255 as pores and grains)
        clip1           = Clip(data, ClipType = 'Scalar', Scalars = ['CELLS', 'im'], Value = 127.5, Invert = 1)
        # 3. make a surface of the remaining grains
        extractSurface1 = ExtractSurface(clip1)
        # 4. and triangulate it for stl export
        triangulate1    = Triangulate(extractSurface1)
    
        # 5. finally save it as an stl file
        SaveData(stl_file, proxy = triangulate1)
    
    # main part
    vti_file = 'porous_model.vti'      # input .vti file
    stl_file = 'porous_model.stl'  # output .stl file
    # call function
    write_stl(vti_file, stl_file)


An easy way to exectute the python scrip is to call it from the python shell in paraview :numref:`fig:paraview_figure_fig`. Within the shell, you first make sure that you are in your case directory; then you make the correct folder to hold the .stl file that we will create (:code:`./constant/triSurface`).

.. code:: python

    import os 
    os.chdir("your_case_directory")
    os.mkdir('constant/triSurface')
    exec(open("your_file_name.py").read())


This will create a .stl file, which we will use in the meshing process.

If you used the paraview GUI, you might have to do some cleanup, so that the files are in the correct location. Move the vti file into :code:`./geometry` and the stl file into :code:`./constant/triSurface`, where openFOAM expects it. This assumes that you are in the case directory.

.. code:: bash

    mv ./porous_model.vti ./geometry/
    mv ./porous_model.stl ./constant/triSurface/

.. admonition:: Doing it the proper way

    The more elegant way would have been to avoid the in-build paraview shell and do everything in a jupyter notebook or a stand-alone python file. Unfortunately, installing the :code:`paraview.simple` module can be a pain - and even the paraview conda package is incompatible with other packages (like vtk, which we will need later in the class). 

    Here is a way to get this to work: make a clean conda environment that has paraview and some other useful things, activate the base environment, startup jupyter and choose the newly created py3_htf_paraview kernel for a new notebook, and finally copy and paste the code above into notebook. Try it!

    .. code-block:: bash

        conda create -n py3_htf_paraview python=3 numpy pandas matplotlib paraview scipy ipykernel

        conda activate base 
        jupyter notebook


OpenFOAM case
-------------

Making the mesh
^^^^^^^^^^^^^^^

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


Notice how we set the vertical and horizontal extents to 1196 and 1494, which is the pixel resolution of the image. We will scale it later to physical dimensions.

Now have a look at :code:`system/snappyHexMeshDict`. We will not go into details here, just explore the general structure yourself if you are interested using the resources linked above.

Time to make the mesh! Run each of the steps below individually and check out the results in paraview.

.. code:: bash

    blockMesh
    snappyHexMesh -overwrite
    checkMesh -allTopology -allGeometry
    transformPoints -scale "(1e-6 1e-6 1e-6)"

The final conversion turns everything in micrometer (:math:`10^{-6} m`). Check out the final mesh in paraview!


Boundary conditions
^^^^^^^^^^^^^^^^^^^

Next we need to set boundary conditions. Our "rock" sample is 1196 and 1494 :math:`\mu m` and we will apply a constant pressure drop of 1 Pa. We will further use openFOAM's **simpleFoam** solver, which resolves incompressible steady-state flow. One important thing to remember is that openFOAM uses the *kinematic* pressure in incompressible simulations, which is the pressure divided by density. If we assume that water with a density of :math:`1000kg/m^3` is flowing through our example, we will need to set a constant pressure value of :math:`P_{left}=0.001 m^2/s^2`. The pressure on the right-hand side will be set to zero. 

Open the file p inside the 0 directory from your local left-hand shell.

.. code-block:: bash 

    code 0/p

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 18,30,31,39,40,44,48
    :name: lst:DRP:bc_p
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
        location    "0";
        object      p;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 2 -2 0 0 0 0];

    internalField   uniform 0;

    boundaryField
    {
        top
        {
            type            symmetryPlane;
        }
        inlet
        {
            type            fixedValue;
            value           uniform 0.001;
        }
        bottom
        {
            type            symmetryPlane;
        }
        outlet
        {
            type            fixedValue;
            value           uniform 0;
        }
        solids
        {
            type            zeroGradient;
        }
        frontAndBack
        {
            type            empty;
        }
    }

The boundary conditions are again set for the patches that were defined in the blockMeshDict. Notice how the sides have constant values. The top and bottom are symmetry planes and the internal boundaries with the solid grains get a :code:`zeroGradient` condition. As we are performing a 2-D simulations, the frontAndBack faces are set to :code:`empty`. 

Remember that units are set by the dimensions keyword. The entries refer to the standard SI units [Kg m s K mol A cd]. By having a 2 in the second and third column, we get the correct units for the kinematic pressure. 

We also need to set boundary conditions for the velocity.

.. code-block:: bash 

    code 0/U

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 30, 34, 42
    :name: lst:DRP:bc_u
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
        class       volVectorField;
        location    "0";
        object      U;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 1 -1 0 0 0 0];

    internalField   uniform (0 0 0);

    boundaryField
    {
        top
        {
            type            symmetryPlane;
        }
        inlet
        {
            type            zeroGradient;
        }
        bottom
        {
            type            symmetryPlane;
        }
        outlet
        {
            type            zeroGradient;
        }
        solids
        {
            type            noSlip;
        }
        frontAndBack
        {
            type            empty;
        }
    }

The top and bottom are again symmetry planes and the :code:`zeroGradient` condition ensure that fluids can freely flow in and out. The contacts with the solid grains get :code:`noSlip`, which makes the velocity go to zero; frontAndBack are again set tp :code:`empty`.


Transport properties
^^^^^^^^^^^^^^^^^^^^

Final step is to set the remaining free parameter (viscosity) in :eq:`eq:mom_con`. 

.. code-block:: bash 

    code constant/transportProperties

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 17-19
    :name: lst:DRP:tp
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
        object      transportProperties;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    transportModel  Newtonian;

    nu              [0 2 -1 0 0 0 0] 1e-06;

    // ************************************************************************* //


Again, in incompressible simulations (where we can divide the momentum balance equation by density) the kinematic viscosity is used, which has units of :math:`m^2/s`. The value of :math:`1e{-6} m^2/s` with an assumed density of :math:`1000kg/m^3` is :math:`1 mPa` - a typical value for water at room temperature. 


Case control
^^^^^^^^^^^^

Finally, we need to set some control parameters in :code:`system/controlDict`. Open it and explore the values.

.. code-block:: bash 

    code system/controlDict

.. code-block:: foam 
    :linenos:
    :emphasize-lines: 16, 37-38
    :name: lst:DRP:cdict
    :caption: controlDict of the DRP case.

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
        object      controlDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    application     simpleFoam;

    startFrom       startTime;

    startTime       0;

    stopAt          endTime;

    endTime         1500;

    deltaT          1;

    writeControl    timeStep;

    writeInterval   1500;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  9;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;

The solver is a stead-state solver, so that the settings are quite simple with respect to transient simulations.


Running the case
^^^^^^^^^^^^^^^^
Now we are finally ready to run the case! Just type this into your docker shell:

.. code-block:: bash 

    ./run.sh

or, if you did the steps above by hand:

.. code-block:: bash 

    simpleFoam


Notice how one new directory is appearing, which contains the steady-state solution. Check that the solution has converged by looking at the :code:`log.simpleFoam` file.

.. code-block:: bash

    tail -15 tail -15 log.simpleFoam

    Time = 537

    smoothSolver:  Solving for Ux, Initial residual = 9.48718063e-10, Final residual = 9.48718063e-10, No Iterations 0
    smoothSolver:  Solving for Uy, Initial residual = 9.93804308e-10, Final residual = 9.93804308e-10, No Iterations 0
    GAMG:  Solving for p, Initial residual = 9.71539045e-09, Final residual = 9.71539045e-09, No Iterations 0
    time step continuity errors : sum local = 7.1290424e-07, global = -1.41776643e-07, cumulative = 3.8954571
    ExecutionTime = 30.24 s  ClockTime = 30 s


    SIMPLE solution converged in 537 iterations

    End 


And explore the solution in Paraview!

.. figure:: /_figures/paraview_DRP_solution.*
   :align: center
   :name: fig:porousModel_solution_fig
   :figwidth: 100%

   Predicted velocity field.


Post-processing / Effective permeability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Time to get back to our initial objective. What is the predicted permeability of our sample? We need to postprocess our solution and solve :eq:`eq:darcy_perm`, which requires are few intermediate steps. First, we need to compute the cell volumes, so that we can integrate the velocity over the total volume. Here comes openFOAM's `postProcess <https://cfd.direct/openfoam/user-guide/v7-post-processing-cli/>`_ functionality handy. Just do this:

.. code-block:: bash 

    postProcess -func writeCellVolumes


This will create a new variable :math:`V` in in the output directories, which is the volume of each cell. Next we can save the full solution as a vtk file, so that we can postprocess it with python.

.. code-block:: bash 

    foamToVTK -useTimeName -latestTime -poly


Check that a new :code:`VTK` directory was created. Now we are ready for some more fancy post-processing using python.

Here is an example script for computing permeability. You can turn it into a jupyter notebook or save it as a .py file into your case directory and run it from there.


.. code-block:: python 

    import vtk
    import numpy as np 
    from vtk.util import numpy_support as VN
    import matplotlib.pyplot as plt 

    # load VTK data
    vtkFile = 'VTK/DRP_permeability_2D_537.vtk'
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtkFile)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()


    # extract velocity arrays 
    Vcells = VN.vtk_to_numpy(data.GetCellData().GetArray('V'))
    U = VN.vtk_to_numpy(data.GetCellData().GetArray('U'))
    Umag = np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
    Ux = U[:,0]
    Uy = U[:,1]
    Uz = U[:,2]

    # calculate volume fluxes
    qx = []
    qy = []
    for c in range(len(Vcells)):
        q1 = Vcells[c] * Ux[c]
        q2 = Vcells[c] * Uy[c]
        qx.append(q1)
        qy.append(q2)
        
    # simulation parameters
    DP = 1                       # pressure drop [Pa]
    nu = 1e-06                   # kinematic viscosity [m²/s²]
    rho = 1000                   # density [kg/m³]
    mu = nu * rho                # dynamic viscosity [kg/(m*s)]
    dx = 0.001196                # model length x [m]
    dy = 0.001494                # model width y [m]
    dz = 1e-6                    # model thickness z [m]
    A = dy * dz
    V = dx * A

    # calculate Darcy velocity
    U_Darcy_x = np.sum(qx)/V

    kxx = mu * dx * U_Darcy_x/DP

    print('Bulk permeability: %5.3e m2' % kxx)

We made it! :math:`2.71e^{-11} m^2` is our predicted effective permeability of our sample. 