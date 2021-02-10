Getting started with OpenFoam
=============================

Why use OpenFoam
----------------
`OpenFoam <http://www.openfoam.org>`_ is an open-source collection of over 200 programs for multi-physics numerical simulations. It's mainly used in the field of Computational Fluid Dynamics. It’s written in C++ and is fully parallelized, treats 3D geometries, is easily modifiable, and has an has active user community.

OpenFoam comes in different flavors and the main forks are:

    * `The Foundation version <https://openfoam.org>`_
    * `The ESI version <https://openfoam.com>`_
    * `The foam-extend development versions <http://sourceforge.net/projects/foam-extend/>`_

The downside is that some of these versions are incompatible with each other and features added to one of the forks may not work in the others. We provided docker image uses the Foundation version.


The OpenFoam workflow
---------------------
A typical workflow involves:

    * **Meshing** of the computational domain either OpenFoam tools like blockMesh or snappyHexMesh or using 3rd party tools
    * Applying **boundary and initial conditions**
    * **Solving** the case using one of OpenFoam's solvers like icoFoam, simpleFoam, etc.
    * **Post-processing** the results, often done in paraview. 

In the next chapter, we will have detailed look at a first example case and solve our first problem.

.. admonition:: Tips for Post-processing 

    Make animation using :code:`ffmpeg` and Paraview.

    1. Install :code:`ffmpeg` in the docker container by running the following command.

    .. code-block::  python
        
        sudo apt-get install ffmpeg

    2. Save images in **the shared folder** using Paraview.

    .. figure:: /_static/video/Paraview_saveanimation.png
        :align: center
        :width: 50 %


    3. Convert images (e.g. *.png) to animation (*.mp4). Run the following command in the docker container, assuming the images are save in the :code:`~/HydrothermalFoam_runs/animation` folder and the filename is :code:`T`.

    .. code-block::  python
        
        cd ~/HydrothermalFoam_runs/animation
        ffmpeg -r 10 -f image2 -s 3002x1526  -i T.%4d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p T.mp4

    Please note the filename pattern :code:`T.%4d.png`, e.g. T.0001.png, T.0098.png, ... Make sure you set a correct pattern according to your own case.