.. include:: /include.rst_

.. _L_FVM_Exercise1:

Exercise 1
===========


.. admonition:: Goals

   - explore FVM theory
   - link FVM theory and OpenFOAM solver implementation
   - step by step exploring details of the solver and solving process


Step 1
------

Create a 2D regular box mesh using :code:`blockMesh` utility (:download:`Regular box case <cases/regularBox.zip>`).

.. figure:: /_figures/mesh_FVM_regularBox.*
   :align: center
   :name: setup_l5_1_fig
   :figwidth: 100%

   2D regular mesh structure and topology information of OpenFOAM polyMesh (See :ref:`/lectures/L_FVM/cases/regularBox/jupyter/VisualizeResults.ipynb#1.-Read-and-plot-mesh-information` in the notebook).



.. toctree::
    :maxdepth: 2

    cases/regularBox/jupyter/VisualizeResults.ipynb


