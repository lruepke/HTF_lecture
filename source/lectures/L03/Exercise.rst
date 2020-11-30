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

Remember how the vertice counting is done! If you don't, check last weeks lecture. 

.. figure:: /_figures/exercise_2layer.*
   :align: center
   :name: fig:exercise_2layer_fig

   Starting point for the two layer system.

The easiest way is start with the blockMeshDict file of the 1-layer case and then add a second layer on top. A good way of getting this right, is to draw the node numbering and then add the vertices to the blockMeshDict file. Notice how we give labels layer_bot and layer_top to the two blocks.

Assign properties
^^^^^^^^^^^^^^^^^

Next we need to set the permeability for the different zones. This can be achieved with :code:`setFields` tool. The :code:`setFields` tool requires a dictionary file :code:`setFieldsDict` that tells it what to do. This file resides in the system folder. Here is a possible listing, just copy it and save it to your system folder.

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

Assign properties
^^^^^^^^^^^^^^^^^
Check all other settings, run the case, and hope for the best!