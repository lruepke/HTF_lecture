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
Again, you can use command of :code:`tree` to explore what changes in the case folder, you can get something like :numref:`lst:flange_mesh:tree`, please notice the new files emphasized in the lines 7-13. 
The :code:`polyMesh` directory contains all the mesh information.

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

    :code:`touch a.foam` means create a empty file named :code:`a.foam`, this could not work in windows system.

.. only:: html

    .. tabs::

        .. tab:: Mesh of flange

            .. figure:: /_figures/flange_mesh.png
                :align: center

                Mesh of the first case - flange.

        .. tab:: Interactive

            .. raw:: html
                
                <iframe src="../../_static/vtk_js/index.html?fileURL=flange_mesh.vtkjs" width="100%" height="500px"></iframe>
        
Boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: /_static/cases/flange/0/T
   :language: c++
   :linenos:
   :lines: 8- 
   :emphasize-lines: 10,12
   :caption: :code:`0/T`.
   :name: lst:flange:0:T

.. admonition:: What do we care about the input file ?

    - boundary/initial conditions are in the 0 directory
    - checkout the 0/T file
    - Flange is initially at 273K
    - patch2 and patch4 have the boundary conditions
    - Units! **[Kg m s K mol A cd]**, unit in the line 10 of :numref:`lst:flange:0:T` is read as :math:`kg^0 \cdot m^0 \cdot s^0 \cdot K^1 \cdot mol^0 \cdot A^0 \cdot cd^0 = K`, it is the unit of temperature.

Transport properties
^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: /_static/cases/flange/constant/transportProperties
   :language: c++
   :linenos:
   :lines: 8- 
   :emphasize-lines: 0
   :caption: :code:`constant/transportProperties`.
   :name: lst:flange:constant:transportProperties

The diffusivity :code:`DT` is specified in dictionary file of :code:`constant/transportProperties`, unit is :math:`kg^0 \cdot m^2 \cdot s^{-1} \cdot K^0 \cdot mol^0 \cdot A^0 \cdot cd^0 = m^2/s`, see :ref:`lecture1_flange_physics`.

Solver and controlDict
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: /_static/cases/flange/system/controlDict
   :language: c++
   :linenos:
   :lines: 8- 
   :emphasize-lines: 11, 13, 15, 17, 19, 21, 25
   :caption: :code:`system/controlDict`.
   :name: lst:flange:system:controlDict

.. admonition:: What should we care about ?

    - solver
    - start/end times
    - time step 
    - output writing 

.. _lecture1_flange_physics:

Physics behind the case 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::
    :label: eq:conti
    
    \frac{\partial T}{\partial t} = \nabla \cdot D_T \nabla T

where :math:`D_T` is the thermal diffusivity.

Results and post-processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. only:: html

    .. tabs::

        .. tab:: Solution in ParaView

            .. figure:: /_figures/flange_solution_paraview.png
                :align: center

                Look at the solution with Paraview
                
        .. tab:: Animation of T evolution

            .. raw:: html

                <video width=100% autoplay loop>
                <source src="../../_static/video/flange_T.mp4" type="video/mp4">
                Your browser does not support HTML video.
                </video>


