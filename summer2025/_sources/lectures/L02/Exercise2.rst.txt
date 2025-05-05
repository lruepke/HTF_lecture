.. include:: /include.rst_

.. _L02_Exercise2:

Excercise 2
===========

After we have mananged to solve for 2-D cavity flow, it's time to try 3-D. And parralel proccessing!

3-D cavity flow
-----------------------------------
First, try to make a 3-D simulations. Start by copying the case to a new folder and to clean the case

.. code-block:: bash 
    
    cp -r cavity2D cavity3D
    ./clean.sh

.. admonition:: Now these things have to happen

    - make the mesh 3-D by having cells in the z-direction
    - make sure the the front and back sides are not set to "empty" anymore
    - set proper boundary conditions for front and back
    - run it and hope for the best!


Getting fancy - let's make a parallel run
-----------------------------------------

Three-dimensional simulations need a lot of computing resources and thankfully OpenFoam is fully parallelized and very efficient in making parallel runs. The standard parallelization technique in OpenFoam is called "domain decompositon". This means that the compuational domain is split into sub-domains and each sub-domain is send to one of the computional cores. As always, all of this is controlled via a dictionary file. This file is called :code:`system/decomposeParDict` .

To get started we make a copy of our 3-D case.

.. code-block:: bash 
    
    cp -r cavity3D cavity3D_par
    ./clean.sh

Unfortunately, theree is no file :code:`system/decomposeParDict` in our example case... Fortunately, there is probably an example file somewhere in all the provided tutorials. A simple way to find and copy one is this:

.. code-block:: bash 
    
    find  $FOAM_TUTORIALS -name "decomposeParDict"
    /opt/OpenFOAM-9/tutorials/heatTransfer/chtMultiRegionFoam/heatedDuct/system/fluid/decomposeParDict
    /opt/OpenFOAM-9/tutorials/heatTransfer/chtMultiRegionFoam/heatedDuct/system/decomposeParDict
    /opt/OpenFOAM-9/tutorials/heatTransfer/chtMultiRegionFoam/heatedDuct/system/metal/decomposeParDict
    /opt/OpenFOAM-9/tutorials/heatTransfer/buoyantSimpleFoam/externalCoupledCavity/system/decomposeParDict
    /opt/OpenFOAM-9/tutorials/heatTransfer/buoyantSimpleFoam/iglooWithFridges/system/decomposeParDict
    cp /opt/OpenFOAM-9/tutorials/heatTransfer/buoyantSimpleFoam/iglooWithFridges/system/decomposeParDict $HOME/HydrothermalFoam_runs/cavity3D_par

Now open this file and edit it; it should look like this

.. code-block:: foam 
    :name: lst:3dcavity_par
    :emphasize-lines: 18, 24
    :linenos:
    :caption: Domain decomposition and parallel processing.  

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
        object      decomposeParDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    numberOfSubdomains 4;

    method          scotch;

    simpleCoeffs
    {
        n               (2 1 2);
        delta           0.001;
    }
    // ************************************************************************* //

The domain is split into four equally sized blocks. Again details on this can be found in the `OpenFoam User Guide <https://cfd.direct/openfoam/user-guide/>`_.

.. figure:: /_figures/decomposed_mesh.*
   :align: center
   :name: fig:decomposed_mesh_fig


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

