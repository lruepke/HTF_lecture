.. include:: /include.rst_

.. _L03_Exercise:

Excercises
==========
In this exercise we will investigate how permeability contrast affect hydrothermal flow and vent temperatures. We will first modify the Regular2DBox case to have constant temperature boundary condition and afterwards add a second layer to our mesh, so that we can explore permeability constrasts.

Constant temperature 1-layer system
-----------------------------------
.. admonition:: What's needed for the constant temperature case?

    - make a copy of the case directory
    - change the temperature boundary condition.

.. figure:: /_figures/RegularBox2D_constantT.*
   :align: center
   :name: fig:RegularBox2D_constantT_fig

   Constant temperature boundary condition.

Constant temperature 2-layer system
------------------------------------
Now we will explore permeability contrasts!

.. admonition:: What's needed for the two-layer case?

    - make a copy of the 1-layer case directory
    - change the mesh to have two layers
    - assign different permeabilities to the two layers


Change the mesh
^^^^^^^^^^^^^^^

Remember how the vertice counting is done! If you don't, check the info box. 

.. figure:: /_figures/exercise_2layer.*
   :align: center
   :name: fig:exercise_2layer_fig

   Starting point for the two layer system.

The easiest way is start with the blockMeshDict file of the 1-layer case and then add a second layer on top. A good way of getting this right, is to draw the node numbering and then add the vertices to the blockMeshDict file. Notice how we give labels layer_bot and layer_top to the two blocks.

.. admonition:: Order of vertices

    .. figure:: /_figures/vertices_order.*
        :align: center
        :name: fig:vertices_order_fig


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


Assign properties
^^^^^^^^^^^^^^^^^

Next we need to set the permeability for the different zones. This can be achieved with :code:`setFields` tool. The :code:`setFields` tool requires a dictionary file :code:`setFieldsDict` that tells it what to do. This file resides in the system folder. Here is a possible listing, just copy it and save it to your system folder.

.. tab:: zoneToCell 

    .. code-block:: foam 
        :linenos:
        :emphasize-lines: 24-34
        :name: lst:2dbox:setfdict
        :caption: Use the setFieldsDict to assign different permeabilities.

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
            location    "system";
            object      setFieldsDict;
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


        defaultFieldValues
        (
            volScalarFieldValue permeability 1e-14
        );

        regions
        (
            zoneToCell
            {
                name "layer_top";
                fieldValues
                (
                    volScalarFieldValue permeability 1e-13
                );
            }
        );

        // ************************************************************************* //

.. tab:: boxToCell 

    .. code-block:: foam 
        :linenos:
        :emphasize-lines: 24-34
        :name: lst:2dbox:boxToCell
        :caption: Use boxToCell to simply set permeability for two-layer model.

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
            location    "system";
            object      setFieldsDict;
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


        defaultFieldValues
        (
            volScalarFieldValue permeability 1e-14
        );

        regions
        (
            boxToCell
            {
                box (0 -500 0) (1000 0 1); //(xmin,ymin,zmin) (xmax,ymax,zmax)
                fieldValues
                (
                    volScalarFieldValue permeability 1e-13
                );
            }
        );

        // ************************************************************************* //



.. tip::

    Alternately, one can also set permeability of two layers by using :code:`boxToCell` based on one layer mesh (see :numref:`lst:2dbox:boxToCell`), which is similar to :code:`zoneToCell` shown in :numref:`lst:2dbox:setfdict`.

Run the case
^^^^^^^^^^^^^^^^^
Now we are using the :code:`setFields` utility to set the permeability. Therefore we need to change the :code:`run.sh` script to also include the setFields command.

.. code-block:: bash 
    :linenos:
    :emphasize-lines: 0
    :name: lst:2dbox:setrun
    :caption: Modified run.sh file that also includes the setFields command.

    #!/bin/sh
    cd ${0%/*} || exit 1    # Run from this directory

    # Source tutorial run functions
    . $WM_PROJECT_DIR/bin/tools/RunFunctions

    application=`getApplication`

    ./clean.sh
    runApplication blockMesh
    cp 0/permeability 0/permeability.orig
    runApplication setFields
    runApplication $application


.. tip::
    Notice that we have also included the statement :code:`cp 0/permeability 0/permeability.orig` into the :code:`run.sh` script. The :code:`setFields` command writes mesh-dependent information into the :code:`permeability` file, which causes problems when we want to change the mesh. Compare the :code:`permeability` file before and after running the :code:`setFields` command. To preserve the old file, we make a copy.

.. only:: html

   Results of the two layer model. The highly permeable upper layer results in cold fluid entrainment and cooling due to mixing.

   .. raw:: html

      <video width=100% autoplay muted controls loop>
         <source src="../../_static/video/layered_movie.mp4" type="video/mp4">
         Your browser does not support HTML video.
      </video>
