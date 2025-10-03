.. include:: /include.rst_

Lecture overview
================

Introduction
----------------

Now that we have OpenFOAM set up, it's time to think about our scientific questions. How can we use numerical simulations to investigate the remote processes within submarine hydrothermal systems. 

Below is a cartoon illustrating the basic processes in a submarine hydrothermal system :cite:`ingebritsen2010numerical`. Cold seawater penetrates the oceanic crust through a network of fractures and porous rocks. As it descends, it heats up due to the geothermal gradient, potentially interacts with a driving magmatic intrusions, and eventually reaches temperatures high enough to leach metals from the surrounding rock. The hot, buoyant fluid then rises back to the seafloor, where it discharges at so-called black smokers. 

.. figure:: /_figures/hts_cartoon.*
   :align: center
   :name: fig:hts_cartoon

   Cartoon illustrating flow in a submarine hydrothermal circulation system :cite:`ingebritsen2010numerical`.


The animation below shows a numerical simulation of such a system using HydrothermalFoam :cite:`hasenclever2014hybrid`. It illustrates possible fluid pathways and temperatures at the fast spreading East Pacific Rise at 9N. The resolved spatial scale is on the order of a few kilometers, while the temporal scale is on the order of several thousand years.

.. only:: html

   A numerical of hydrothermal flow at the East Pacific Rise at 9N by :cite:`hasenclever2014hybrid`.

   .. raw:: html

      <video width=100% autoplay muted controls loop>
         <source src="../../_static/video/Hasenclever_EPR_small.mp4" type="video/mp4">
         Your browser does not support HTML video.
      </video>


In this course, we will learn how to set up and run such simulations using OpenFOAM and HydrothermalFoam. 