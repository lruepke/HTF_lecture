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

Mesh generation
^^^^^^^^^^^^^^^^^^^^^^^

The mesh of flange case is prepared as Ansys_ format (see line 8 of :numref:`lst:flange:tree`), you can use OpenFOAM mesh convert utility of :code:`ansysToFoam` to convert the ansys format to OpenFOAM format, the command is shown below,

.. code-block:: bash 

    ansysToFoam flange.ans -scale 0.001

The :code:`-scale` option is used to set scale factor of the whole mesh. 
After converting the mesh, you can use Paraview_ to visualize the mesh,

.. code-block:: bash

    touch a.foam
    paraview a.foam 

.. note:: 

    :code:`touch a.foam` means create a empty file named :code:`a.foam`, this could not work in windows system.

.. only:: html

    .. tabs::

        .. tab:: Mesh of flange

            .. figure:: /_figures/flange_mesh.png
                :align: center

                Mesh of the first case - flange.

        .. tab:: Interactive

            Mesh and boundary condions of the 3D pipe model.

            .. raw:: html
                
                <iframe src="../../_static/vtk_js/index.html?fileURL=pipe_3D_mesh.vtkjs" width="100%" height="500px"></iframe>
        
