.. include:: /include.rst_

.. _L05_Exercise1:

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

   - modify the mesh (blockMeshDict)
   - check that gravity is zero (constant/g)
   - update boundary conditions for pressure and temperature (and permeability)
   - use :code:`setFields` utility and system/setFieldsDict to modify the permeability structure
   - the exact of dimensions of the permeability blocks don't matter


.. tip::

    Remember the two-layered convection exercise from lecture 3? There we used the boxToCell function as an alternative to making two mesh zones. Now we can do the same here, just that we need several boxes.

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
