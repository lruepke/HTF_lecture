.. include:: /include.rst_

.. _sec_Exercises:

Exercises
=========

What to do
----------

A first test case - heat diffusion
-----------------------------------------

.. admonition:: What do we need?

    - case files
    - a mesh
    - boundary/initial conditions
    - a solver
    - simulations controls (time step, run time, etc.)
    - post-processing

Prepare case files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first case is the **flange** case from OpenFOAM tutorials, just copy the whole case to your working directory, this can be done using the following command in your terminal (**should be the docker container terminal if you use OpenFOAM in the docker container**).
It should be noted that you have to change directory to your working directory first.

.. code-block:: bash 
    :name: lst:cpFlangeToWorkDir

    cp -rf $FOAM_TUTORIALS/basic/laplacianFoam/flange . 

The tree structure of flange case is shown in :numref:`lst:flange:tree`.

.. code-block:: bash 
    :linenos:
    :emphasize-lines: 0
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

Mesh generation
^^^^^^^^^^^^^^^^^^^^^^^

The mesh of flange case is prepared as Ansys_ format, you can use OpenFOAM mesh convert utility of :code:`ansysToFoam` to convert the ansys format to OpenFOAM format, the command is shown below,

.. code-block:: bash 

    ansysToFoam flange.ans -scale 0.001
