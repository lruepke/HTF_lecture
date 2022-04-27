.. include:: /include.rst_

.. _L03_Exercise:

Exercises
==========
Let's explore the solution in a bit more detail.

Vertical versus horizontal permeability
---------------------------------------
Make a copy of your case and repeat the analysis in the vertical direction. You will need to

    * Change the boundary conditions
    * Modify the post-processing script


Changes to the geometry
------------------------------------
Check if/how the solution changes if you distort the mesh. For example, squeeze it a bit in the vertical direction by modifying the :code:`transformPoints -scale "(1e-6 1e-6 1e-6)"` statement. Also check if the absolute dimensions matter!

