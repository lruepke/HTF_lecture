.. include:: /include.rst_

Lecture overview
================

In this lecture we explore an example problem from "Digital Rock Physics". We will use direct fluiddynamic simulations to assess the hydro-geological properties of porous rocks. As an example problem, we will calculate the effective permeability of a synthetic 2-D representation of the pore space. The steps involved include:

   * Using an image as input and turning it into a pore-space model
   * Meshing the pore space using `openFOAM's snappyHexMesh tool <https://cfd.direct/openfoam/user-guide/v7-snappyHexMesh/#x27-1970005.4>`_  
   * Performing a flow simulations using constant pressure boundary conditions
   * Computing the effective permeability of the sample by evaluating the total flow


