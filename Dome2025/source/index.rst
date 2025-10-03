.. HTF_lecture documentation master file

.. include:: /include.rst_

Introduction to OpenFOAM and HydrothermalFoam
=============================================
Welcome to the DOME short course on hydrothermal flow modeling. The course intends to provide you with a hands-on experience on how we can use numerical methods to investigate the inner workings of submarine hydrothermal systems. For this we will use the popular Computational Fluid Dynamics (CFD) package OpenFOAM_ , which mainly adresses Navier-Stokes-type fluiddynamic problems. In addition, we will have a look at HydrothermalFoam_, a custom solver for flow in porous media developed at GEOMAR, which allows resolving hydrothermal convection in submarine settings.
     

Course content
--------------
These are the main topics:

1. Getting started with OpenFOAM and HydrothermalFoam
2. Direct simulations of flow on the pore scale
3. Upscaled porous flow modeling and submarine hydrothermal systems 


.. admonition:: Open access!

    Note that OpenFOAM_ and HydrothermalFoam_ are open source and can be freely used for your research.


Instructors
-----------

**Prof. Lars Ruepke**

- Email: lruepke@geomar.de
- `Institute webpage <https://www.geomar.de/en/research/fb4/fb4-muhs/research-topics/modelings>`_


Course website
---------------

You can find additional lectures and information on our teaching website:

- Public site : https://lruepke.github.io/HTF_lecture/


Further readings
----------------

There are many good online resources on OpenFOAM, HydrothermalFoam, and numerical modeling in marine geosciences:

- Cyprien Soulaine's teaching material `<http://cypriensoulaine.com>`_
- Official OpenFoam documentation  `<https://cfd.direct/openfoam/documentation/>`_
- Material shared by Tobias Holzmann  `<https://holzmann-cfd.com>`_ 
- OpenFoam code documentation `<https://cpp.openfoam.org/v10/>`_
- CFD Forum `<https://www.cfd-online.com/Forums/>`_


.. toctree::
    :maxdepth: 2
    :caption: Lecture 1

    lectures/L01/Overview   
    lectures/L01/Installation
    lectures/L01/OpenFoam
    lectures/L01/FirstCase

.. toctree::
    :maxdepth: 2
    :caption: Lecture 2

    lectures/L02/Overview
    lectures/L02/Introduction
    lectures/L02/Flow_pore_scale
    lectures/L02/Exercise

.. toctree::
    :maxdepth: 2
    :caption: Lecture 3

    lectures/L04/Overview
    lectures/L04/intro  
    lectures/L04/FirstCase  
    lectures/L04/Exercise  


.. toctree::
   :maxdepth: 2
   :caption: References
   
   refs


    
