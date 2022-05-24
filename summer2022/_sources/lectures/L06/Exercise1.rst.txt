.. include:: /include.rst_

.. _L06_Exercise1:

Exercise 1
===========

We will first explore a simple setup of 2-D flow along a pressure gradient in the presence of permeability contrasts.

.. admonition:: Goals

   - explore porous flow along prescribed pressure gradient
   - use :code:`setFields` utility to modify the permeability structure
   - Understand how permeability contrasts affect flow

Step 1
------

We will explore the following case: A 2-D box of dimensions 1500 m by 1000 m. Flow is driven by a prescribed pressure gradient. Within the box, variations in permeability affect the flow field.

.. figure:: /_figures/lecture_plot.*
   :align: center
   :name: setup_l5_1_fig
   :figwidth: 80%

   Setup of our 2D test case for flow in a porous medium with variable permeability.


A good starting point for setting up new case is always to take an existing similar case. Here we can use one of |foam| benchmark cases. It's located in :code:`$HOME/hydrothermalfoam-master/benchmarks/HydrothermalFoam/1d/h1`.

Copy the case into your working directory.

.. code-block:: bash

    cd $HOME/HydrothermalFoam_runs
    cp -r /home/openfoam/hydrothermalfoam-master/benchmarks/HydrothermalFoam/1d/h1 ./h1_class

Step 2
------

Now we need to modify the case. These are the necessary steps:

.. admonition:: Goals

   - modify the mesh by chnanging :code:`blockMeshDict`
   - check that gravity is zero in :code:`constant/g`
   - update boundary conditions for pressure and temperature (and permeability)
   - use the :code:`setFields` utility and the controlling :code:`system/setFieldsDict` to modify the permeability structure
   - the exact of dimensions of the permeability blocks don't matter
   - modify :code:`system/controlDict` so that we only make one time step. Easiest way is to set endTime, deltaT and writeInterval all to 864000.


.. tip::

    You might have noticed that the :code:`benchmarks/HydrothermalFoam/1d/h1` case does not contain a :code:`setFieldsDict`. This is common "problem" that a dictionary file is missing when starting from an existing case file. You can easily fix this by copying the respective dictionary from another case and then modify it. If you don't know which case has the missing dictionary, you can just search for it (inside the docker container). For example like this:

    .. code-block:: bash

        find $HOME -type f -name "setFieldsDict"
 
    Now that we have the setFieldsDict, we need to modify it. Do you remember :ref:`L03_properties` part the two-layered convection :ref:`L03_Exercise` in lecture 3? There we used the boxToCell function within :code:`setFieldsDict` as an alternative way of making a layerd permeability structure. Now we can do the same here, just that we need several boxes.

    This is how the code snippit looked like in lecture 3; take it as a starting point!

    .. code-block:: bash

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

    You can also play with different shapes, like cylinder or rotated box. Check the documentation!

Step 3
------

Explore the results in paraview!

.. figure:: /_figures/hetero_perm_flow.*
   :align: center
   :name: flow_l5_1_fig
   :figwidth: 80%

   Porous flow in the presence of permeability variations.

