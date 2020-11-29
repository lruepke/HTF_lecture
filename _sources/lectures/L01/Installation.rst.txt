.. include:: /include.rst_

Installation guide
==================

We will be using the Foundation version of OpenFOAM_ and our own hydrothermal flow solver HydrothermalFoam_ , which is based on OpenFoam. Detailed installation instructions can be found on the webpages of HydrothermalFoam. During this class we will rely on a virtualized environment and use a pre-configured docker container. If you don't like virtual environments, you will need to go through the "install from source" instructions of HydrothermalFoam and read the installation instructions of OpenFoam. 

.. _sec:install_docker:

Docker
------
We will rely on docker_ and use a virtual Ubuntu Linux installation that has OpenFoam and HydrothermalFoam installed. For this purpose we have modified the official OpenFoam docker image, so that it also contains all components of HydrothermalFoam. The way docker works is that a so-called docker image is downloaded and from this image a so-called docker container is built, which is the actual virtual environment. These are the necessary steps:

1. **Install** `docker desktop <https://www.docker.com/products/docker-desktop>`_ and keep it running.
2. **Register** an account with `docker hub <https://hub.docker.com/>`_, a marketplace for docker images.
3. **Log-into** your newly created account from your docker desktop application.

Now we are ready to pull the provided docker image, which can be found here `zguo/hydrothermalfoam <https://hub.docker.com/repository/docker/zguo/hydrothermalfoam>`_.

Open a shell (powershell under windows or a terminal under MacOS) and pull the image by typing this command.

.. code-block:: bash

      docker pull zguo/hydrothermalfoam

This can take a while. After the download is finished, you can build the docker container. Use this command:

.. code-block:: bash

    docker run -it -d --name hydrothermalfoam --workdir="/home/openfoam" - v="$HOME/HydrothermalFoam_runs":"/home/openfoam/HydrothermalFoam_runs" zguo/hydrothermalfoam

This complicated looking command builds the docker container and creates a “shared” folder called HydrothermalFoam_runs in your home directory (this is what the -v option does). If you want it somewhere else, you can change this "$HOME/HydrothermalFoam_runs" part of the above statement. Please do not uses path and directory names with spaces or special characters in it (if you are on windows).

.. tip::
    Just in case you do not like to type statement that you don’t fully understand – this is what the other options do:
        * -it allocate a pseudo-TTY connected to the container’s stdin; creating an interactive bash shell in the container (plain language, this allows you to interact with (type inside) the container)
        * -- name hydrothermalfoam give the container its name
        * --workdir is the path to the working/home directory inside the container (don’t change it!)
        * -v makes a shared file system
        * -d run container in background (detached)


.. _sec:install_vscode:

Visual Studio Code
------------------
Working with OpenFoam involves a lot of editing of text files and you can use your favorite text editor for this. However, we strongly recommend to use `Microsoft's Visual Studio Code <https://code.visualstudio.com/>`_. 

Open Visual Studio Code and open two shells side-by-side at the bottom (powershell on Windows and terminal on Mac). We will use the convention that we use the terminal to the right-hand side as our “virtual machine (docker) shell” and the terminal to the left as our local shell.

Start and attach the docker container to the right-hand-side shell by typing

.. code-block:: bash

    docker start hydrothermalfoam
    docker attach hydrothermalfoam

If you want things to be pretty, try starting a z-shell (instead of the default bash):

.. code-block:: bash

    zsh

You can stop the container whenever you want and your files will remain intact, just like in any other virtual machine. To stop the container type

.. code-block:: bash
    
    docker stop hydrothermalfoam

.. tip::
    It really helps to have syntax highlighting when manipulating OpenFoam files. To get it search for openfoam in the Extensions Tab of Visual Studio Code and you will find an extension provided by Zhikui Guo. Install the plug-in and copy the file associations into the settings.json file. To do so, type CMD/CTRL + SHIFT + P and search for "Open Settings (JSON)".

.. _sec:install_paraview:

Paraview
--------
We will use Kiteware's `Paraview <paraview.org>`_ to visualize the modeling results. Download and install the latest version. It also really helps to add paraview to the system path, so that it can be called from the command line.