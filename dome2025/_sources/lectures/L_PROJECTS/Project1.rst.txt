.. include:: /include.rst_

.. _L_PROJECTS_Project1:

Project 1: Flow along a detachment fault
===========================================

.. admonition:: Project objectives

    - Explore hydrothermal flow along fault zones
    - Find which faults capture a rising plume and which do not
    - Understand how fault zone properties affect vent temperature


The first project builds upon lecture 3 on :ref:`fault_controlled_systems`. We will explore under which conditions hydrothermal upflow can be "captured" by an inclined fault zone of enhanced permeability. 

.. figure:: /_figures/andersen_fault_2017.*
   :align: center
   :name: fig:fault_2017_fig
   :figwidth: 70%

   Figure taken from :cite:`Andersen2017` ; preferential fluid flow along a fault zone.



Deliverables
-----------------------

The goal of this project is to provide work on these topics

1. **Plume captureing**, provide a sequence of 2D simulations that illustrate what can "happen" for different fault widths and fault permeabilities (check :numref:`andersen_fig`)
2. **Systematics**, provide a "phase diagram" for fault width and fault permeability that illustrates the regime when the pume is captured by the fault zone.
3. **Vent temperature**, provide information on the relationship between the fault's ability to capture a plume and the resulting vent temperature.

The results should be delivered in a short powerpoint presentation.

Starting point
--------------

Case file
^^^^^^^^^^
Download the basic case file from (:download:`Fault Flow Model <cases/fault_flow.zip>`).. It includes a basic setup for a fault-controlled hydrothermal system.

Meshing software
^^^^^^^^^^^^^^^^

The example case uses a meshing software called Gmsh_ to create a mesh that resolves the fault zone. Gmsh is already installed within your docker containers, so there should be no need for installing additional software. The gmsh input file is in :code:`gmsh/make_mesh.geo`. 

.. figure:: /_figures/gmsh_screen.*
   :align: center
   :name: fig:gmsh_screen_fig

   Mesh created with gmsh.

Have a look at :code:`gmsh/make_mesh.geo` and the Gmsh_ . Check that you (more or less) understand what's happening and try to change things like fault width and and angle.

Post-processing
^^^^^^^^^^^^^^^

There are many ways of extracting information from a completed run. For this project, you will want to know how high the vent temperature is and where the venting occurs, so that you can evaluate if the plume was "captured" by the fault. You can either do this in paraview by using, for example, the "Plot Over Line" filter in Paraview. An alternative is to use the built-in postprocessing function of OpenFoam. 

.. code::

    postProcess -func sampleDict -latestTime

It will read the :code:`system/sampleDict` file (have a look!) and extract T along a line.

More powerful is, of course, to use python, like in :ref:`L04_PostProcess`. 


Background reading
------------------

    * :cite:`andersen2015` illustrates the problem for the Logatchev field
    * :cite:`driesner2010` explains relationship between vent temperature and permeability
    * :cite:`deMartin2007` background paper on detachment faulting at the TAG vent field
    * :cite:`McCaig2007` :cite:`McCaig2010` more paper on fluid flow along detachments
    * :cite:`Escartin2008` classic paper on oceanic detachments

