.. include:: /include.rst_

.. _L02_Exercise2:

Exercise 2
===========

After we have managed to solve for 2-D cavity flow, it's time to try 3-D. And parallel processing!

3-D cavity flow
-----------------------------------
First, try to make a 3-D simulation. Start by copying the case to a new folder and clean it

.. code-block:: bash 
    
    cd $HOME/HydrothermalFoam_runs
    cp -r cavity2D cavity3D
    cd cavity3D
    ./clean.sh

.. admonition:: Now these things have to happen

    - make the mesh 3-D by having cells in the z-direction (:code:`system/blockMeshDict`)
    - make the domain a cube: the cavity is only 0.01 m thick so far. Which vertex coordinates do you have to change?
    - the front and back sides must not be "empty" anymore. There are three places where you have to change this: the patch type in :code:`system/blockMeshDict`, and the boundary conditions in :code:`0/U` and :code:`0/p`
    - set proper boundary conditions for front and back: they are walls now, just like the fixed walls
    - check the Courant number in the log file
    - run it and hope for the best!

.. tip::
    If something goes wrong, read the error message carefully, it usually names the file and the patch with the problem. A frequent error is a patch that is :code:`empty` in one file but not in the others.


Getting fancy - let's make a parallel run
-----------------------------------------

Three-dimensional simulations need a lot of computing resources and thankfully OpenFOAM is fully parallelized and very efficient in making parallel runs. The standard parallelization technique in OpenFOAM is called "domain decomposition". This means that the computational domain is split into sub-domains and each sub-domain is sent to one of the computational cores. As always, all of this is controlled via a dictionary file. This file is called :code:`system/decomposeParDict` .

To get started we make a copy of our 3-D case.

.. code-block:: bash 
    
    cd $HOME/HydrothermalFoam_runs
    cp -r cavity3D cavity3D_par
    cd cavity3D_par
    ./clean.sh

Unfortunately, there is no file :code:`system/decomposeParDict` in our example case... Fortunately, there is probably an example file somewhere in all the provided tutorials. A simple way to find and copy one is this:

.. code-block:: bash 
    
    find  $FOAM_TUTORIALS -name "decomposeParDict"
    ...
    /opt/OpenFOAM-10/tutorials/heatTransfer/chtMultiRegionFoam/heatedDuct/system/decomposeParDict
    /opt/OpenFOAM-10/tutorials/heatTransfer/buoyantFoam/externalCoupledCavity/system/decomposeParDict
    /opt/OpenFOAM-10/tutorials/heatTransfer/buoyantFoam/iglooWithFridges/system/decomposeParDict
    ...
    cp /opt/OpenFOAM-10/tutorials/heatTransfer/buoyantFoam/iglooWithFridges/system/decomposeParDict $HOME/HydrothermalFoam_runs/cavity3D_par/system

Now open this file and edit it. The tutorial file is set up for two subdomains; change it to four subdomains and to the method :code:`simple`, so that it looks like this (the other coefficient blocks in the file are not used and you can delete them):

.. code-block:: foam
    :name: lst:3dcavity_par
    :emphasize-lines: 16, 18, 22
    :linenos:
    :caption: Domain decomposition and parallel processing.

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
        object      decomposeParDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    numberOfSubdomains 4;

    method          simple;

    simpleCoeffs
    {
        n               (2 2 1);
    }

    // ************************************************************************* //

With the method :code:`simple`, the domain is split into equally sized blocks: :code:`n (2 2 1)` means two blocks in the x direction, two in the y direction, and one in the z direction. The product of the three numbers must equal :code:`numberOfSubdomains`. Other methods, e.g. :code:`scotch`, split the domain automatically and do not use :code:`simpleCoeffs`. Details can be found in the `OpenFOAM User Guide <https://doc.cfd.direct/openfoam/user-guide-v10/running-applications-parallel>`_.

.. figure:: /_figures/decomposed_mesh.*
   :align: center
   :name: fig:decomposed_mesh_fig

   A 3-D cavity mesh decomposed into four subdomains. Each colour is the part of the mesh that one processor works on.


Next you need to modify the :code:`run.sh` script to invoke a parallel run. 

.. code-block:: bash 
    
    cp run.sh run_par.sh

It should look like this:

 .. code-block:: bash 
    :linenos:
    :name: lst:2dcav_run:tree
    :caption: The run_par.sh file.

    #!/bin/sh
    cd ${0%/*} || exit 1    # Run from this directory

    # Source tutorial run functions
    . $WM_PROJECT_DIR/bin/tools/RunFunctions

    application=`getApplication`

    ./clean.sh
    runApplication blockMesh
    runApplication decomposePar
    runParallel $application
    reconstructPar

Now it's time to run the case. Results should look like this:


.. figure:: /_figures/3dcavity.*
   :align: center
   :name: fig:3dcavity_fig
   
   Results of the 3-D parallel cavity flow simulation.

