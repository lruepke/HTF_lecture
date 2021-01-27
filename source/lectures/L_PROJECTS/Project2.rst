.. include:: /include.rst_

.. _L_PROJECTS_Project2:

Project 2: How long can a magmatic intrusion sustain high temperature venting?
==============================================================================

. admonition:: Project objectives

    - Explore how long magmatic intrusions can sustain a hydrothermal circulation system
    - Constrain the break-through time of the hydrothermal plume
    - Interpret the results in the context of mineral systems

The second exercise explores how long a magmatic intrusion can sustain a high-temperature hydrothermal system and how much of a time lag exists between the magmatic event and the formation of a vent site. These king of insights are highly relevant to large number of research questions.

Lifetimes of hydrothermal systems and associated mineral systems

For example, sill intrusions into sedimentary basins are thought to be related to the release of large volumes of greenhouse and toxic gases into the atmosphere. 


.. figure:: /_figures/3d_sill.*
   :align: center
   :name: fig:3d_sill_fig

   Figure taken from :cite:`Galerne2019`on sill intrusions and hydrothermal venting in sedimentary basins.


something on hydrothermal systems and sealeve...

Technicalities
----------------

Solver changes
^^^^^^^^^^^^^^^

The cooling of magmatic intrusions is associated with a latent effect of the crystallization process. A convenient way of accounting for this effect is to use an effective specific heat for the solid phase in the energy equation (see e.g. :cite:`Aarnes2010`). 

We can do this by modifying the transient term of the energy equation:

.. math::
    :label: eq:e_eqn_eff
    
    (\varepsilon \rho_f C_{pf} + (1-\varepsilon)\rho_r C_{eff})\frac{\partial T}{\partial t}  = ...

.. math::
    :label: eq:e_eqn_cp1
    
    C_{eff}  = C_p \left( 1+ \frac{L_c}{(T_L -T_S) C_p)}\right) \text{ for } (T_S < T < T_L)

.. math::
    :label: eq:e_eqn_cp2
    
    C_{eff}  = C_P \text{ for } (T < T_S)   

These changes are implemented in a modified solver of HydrothermalFoam, which you can download from here (:download:`HydrothermalSinglePhaseDarcyFoam_Cpr <cases/HydrothermalSinglePhaseDarcyFoam_Cpr.zip>`).

Check the :code:`updateProps.H` and :code:`createFields.H` files, how these changes are implemented.

Setup fields
^^^^^^^^^^^^^^^

In addition to the solver changes, we have to take extra care of the initial conditions. The temperature field is not uniform anymore but we need to set the initial intrusion temperature to a high value according to its cellzone and apply a background temperature gradient (if we think that's useful).

These changes are implemented as codestream statements in :code:`0/T` . Have a look!


Deliverables
-----------------------

The goal of this project is to provide work on these topics

1. **Intrusion driven circulation**, provide a sequence of 2D simulations that illustrate what can "happen" for different intrusion geometries and host rock permeabilities.
2. **Systematics**, provide "phase diagrams" on breakthrough times, the time to vent, and vent durations, the time high-T venting is sustained, for different model parameters like intrusion volume, depth, permeability..
3. **Vent temperature**, provide information on the relationship between the host rock permeability and maximum vent temperature. 

The results showed be delivered in a short powerpoint presentation.


Starting point
--------------

Case file
^^^^^^^^^^
Download the basic case file from (:download:`Intrusion Flow Model <cases/cooling_intrusion.zip>`).. It includes a basic setup for hydrothermal system drive.

Meshing software
^^^^^^^^^^^^^^^^

The example case uses a meshing software called `gmsh  <https://gmsh.info>`_ to create a mesh that resolves the fault zone. Gmsh is already installed within your docker containers, so there should be no need for installing additional software. The gmsh input file is in :code:`gmsh/make_mesh.geo`. 

.. figure:: /_figures/sill_mesh.*
   :align: center
   :name: fig:gmsh_sill_fig

   Mesh created with gmsh.

Have a look at :code:`gmsh/make_mesh.geo` and the `gmsh documentation  <https://gmsh.info>`_ . Check that you (more or less) understand what's happening and try to change things like fault width and and angle.

Post-processing
^^^^^^^^^^^^^^^

There are many ways of extracting information from a completed run. For this project, you will want to know how high the vent temperature is and where the venting occurs, so that you can evaluate if the plume was "captured" by the fault. You can either do this in paraview by using, for example, the "Plot Over Line" filter in Paraview. An alternative is to use the built-in postprocessing function of OpenFoam. 

Also check the :code:`system/controlDict.orig` file; there is codestream section at the bottom that writes out maximum vent temperature to a file names :code:`ventT.txt` . Modify according to your needs.


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

