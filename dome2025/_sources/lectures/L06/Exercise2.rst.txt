.. include:: /include.rst_

.. _L06_Exercise2:

Exercise 2
===========

Ok, time to put our new knowledge to work! Let's explore a more realistic case and investigate how a more permeable fault zone can deviate hydrothermal upflow. You can find background on such a setup in :cite:`andersen2015`.

.. admonition:: Steps

   - check all settings and modify model run time so that the plume arrives at the surface
   - modify (or create) a :code:`setFieldsDict` to create the permeability structure. 
   - use a function called :code:`rotatedBoxToCell` to set the permeability structure
   - run the case and explore how the results change for different fault widths and permeability contrasts.


Step 1
-------

Copy the :code:`$HOME/hydrothermalfoam-master/cookbooks/2d/Regular2DBox/` into your working directory.

.. code-block:: bash

    cd $HOME/HydrothermalFoam_runs
    cp -r /home/openfoam/hydrothermalfoam-master/cookbooks/2d/Regular2DBox/ ./Fault_zone_2D


Step 2
------

Now modify the case and setup the permeability structure. It should look like this:

.. figure:: /_figures/perm_fault.*
   :align: center
   :name: setup_l5_2_fig
   :figwidth: 80%

   Permeability structure of our detachment fault example cause.


We can use the :code:`rotatedBoxToCell` function in the :code:`setFieldsDict` for this. Check the `openfoam documentation <https://cpp.openfoam.org/v11/classFoam_1_1rotatedBoxToCell.html>`_ !

.. tip::

    It's all about rotating the coordinate system to figure out the origin and i,j,k vectors for :code:`rotatedBoxToCell`!

    .. figure:: /_figures/coord_sys.*
        :align: center
        :name: coord_sys_fig
        :figwidth: 40%


You can use the :code:`#calc` macro that we already used in the very first lecture to do the rotation. The :code:`setFieldsDict` should contain line like this:


.. code-block:: foam

    // define variables
    // origin of rotated box
    X0 750;
    Y0 -2800;
    Z0 0;

    width  100; //width of fault zone 
    height 800; // distance to surface of box
    angle 60;  // degrees, counter-clockwise angle between x-axis and fault

    // Example using #calc macro
    angleRadians  #calc "degToRad($angle)";  // Convert angle to radians

    // input for rotatedBoxToCell
    i0 ???;
    i1 ???;
    i2 ???;

    j0 ???;
    j1 ???;
    j2 ???;

    regions
    (
        rotatedBoxToCell
        {
            origin   ( ? ? ?);
            i        ( ? ? ?);
            j        ( ? ? ?);
            k        ( ? ? ?);

            fieldValues
            (
                volScalarFieldValue permeability 1e-13
            );
        }
    );



Step 3
------

Explore the results in paraview! And investigate how the results change for different fault widths and permeability contrasts.

.. figure:: /_figures/flow_fault.*
   :align: center
   :name: flow_fault_fig
   :figwidth: 80%

   Hydrothermal flow along a preferential pathway.



.. tip::

    - Run many simulations for different fault widths and permeabilities
    - Explore how vent temperature changes
    - ...and when the plume is captured by the fault!