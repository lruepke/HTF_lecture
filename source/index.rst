.. HTF_lecture documentation master file, created by
   sphinx-quickstart on Sat Nov 28 21:47:26 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. highlight:: console

Introduction to OpenFOAM and HydrothermalFoam
=============================================

Welcome to the website for the course Numerical Fluid Dynamics I (flow in porous media) taught at Kiel University within the Geophysics master.  The course introduces methods and tools in geophysical fluid dynamics. The goal is to learn about the complete workflow from formulating a geoscientific hypothesis to testing it using numerical techniques. For this we will use the popular Computational Fluid Dynamics (CFD) package `OpenFOAM <http://openfoam.org>`_. In terms of numerical methods, we will learn about finite-differences and finite-volumes including their respective strengths and limitations. Topic-wise, we will focus on porous flow processes with a special emphasis on hydrothermal flow problems like submarine black smoker systems and hydrothermal cooling of magmatic intrusions. For these topics we will use the 3-D hydrothermal flow model `HydrothermalFoam <http://hydrothermalfoam.info>`_, which is also based on OpenFOAM.


Course content
--------------
These are the main topics:

1. Introduction to computational fluids dynamics /  Navier-Stokes, Stokes, Darcy
2. Numerical methods / finite differences and finite volumes (MATLAB/Python)
3. Computational Fluid Dynamics with OpenFOAM 
4. Visualization with Paraview
5. Hydrothermal systems
6. Project work on transport in porous media

Course goals
--------------

1. Develop the ability to perform independent modeling work on geoscientific problems.
2. Learn numerical techniques and of how to solve partial differential equations using numerical methods.
3. Obtain in-depth knowledge of OpenFOAM (problem setup, meshing, solvers, 2D -> 3D, visualization). 
4. Perform independent project work.


Course format
-------------

The majority of this course will be spent in front of a computer working on exercises related to (porous) flow problems in marine geodynamics.

.. code-block:: bash 
    :name: lst:createCaseDir

    mkdir test
    cd !$  



.. admonition:: Open access!

    Note that the materials for this course are **open to everyone**; the course is, however, taught as an on-site class for registered students at Kiel University. 

.. toctree::
    :maxdepth: 2
    :caption: Course information

    general-info/course-details

.. toctree::
    :maxdepth: 2
    :caption: Lecture 1

    lectures/L01/intro
    lectures/L01/exercise
