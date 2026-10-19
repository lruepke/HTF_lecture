.. include:: /include.rst_

.. _Installation guide:

Installation guide
==================

We will be using the Foundation version of OpenFOAM_ (version 10) and our own hydrothermal flow solver HydrothermalFoam_, which is based on OpenFOAM. Everything you need for the class comes in one pre-configured `Docker <https://www.docker.com>`_ container:

* OpenFOAM 10, built from source on Ubuntu 22.04 Linux
* HydrothermalFoam with its solver :code:`HTFoam`, and xThermo, our implementation of the :cite:`Driesner2007` equation-of-state for the system H\ :sub:`2`\ O-NaCl
* Python 3.12 with all packages we use in class (numpy, scipy, matplotlib, pandas, vtk, meshio, porespy, iapws, ...)
* JupyterLab and Visual Studio Code, which you use in your web browser

You only install two programs on your computer: **Docker Desktop**, which runs the container, and **ParaView**, which we use to look at the results.

.. tip::
    Docker images are templates that contain an operating system with installed software. A container is a running instance of an image, a bit like a virtual Linux computer inside your computer. You can stop and restart it, and your files stay where they are.

.. _sec:install_docker:

Docker Desktop
--------------

1. **Install** `Docker Desktop <https://www.docker.com/products/docker-desktop>`_ and start it. It needs to be running whenever you work on the class.
2. **Optional:** register an account at `Docker Hub <https://hub.docker.com/>`_ and log in from Docker Desktop. You can download the course image also without an account.

.. tip::
    The course image needs about 8 GB of disk space. On Windows, Docker Desktop uses WSL 2; just accept the defaults during the installation.

.. _sec:install_paraview:

ParaView
--------
We will use Kitware's ParaView_ to visualize the modeling results. Download and install the latest version for your computer. ParaView runs on your computer, not in the container, and opens the result files from your course folder.

.. _sec:install_container:

The course container
--------------------

1. Create a **course folder** on your computer, for example :code:`HTF_course` in your home directory. Please do not use spaces or special characters in the folder name and in the path to it.
2. Download the file :download:`docker-compose.yml <docker-compose.yml>` and save it into your course folder. It describes the course container.
3. Open a terminal (Terminal on macOS and Linux, PowerShell on Windows), go to your course folder, and start the container:

   .. code-block:: bash

       cd HTF_course
       docker compose up -d

   The first time, this downloads the course image, which takes a while. Afterwards, the container starts within seconds.

4. Open the two tools in your web browser:

   * **Visual Studio Code**: http://localhost:8080
   * **JupyterLab**: http://localhost:8888

Docker created a folder :code:`HydrothermalFoam_runs` in your course folder. Inside the container, the same folder is :code:`/home/openfoam/HydrothermalFoam_runs`. This is where all your work goes: everything you save there is directly visible on your computer (e.g. for ParaView), and everything you put there on your computer is visible in the container.

.. warning::
    Only the files in :code:`HydrothermalFoam_runs` are stored on your computer. Everything else in the container is lost when the container is removed (e.g. by :code:`docker compose down` or when the course image is updated). So always work in :code:`HydrothermalFoam_runs`.

Daily use
^^^^^^^^^
Run these commands in a terminal in your course folder:

.. code-block:: bash

    docker compose up -d      # start the container
    docker compose stop       # stop the container when you are done

You can also start and stop the container in Docker Desktop.

.. figure:: /_figures/docker_desktop.png
    :align: center
    :figwidth: 70%

    Docker Desktop shows your images and containers.

.. _sec:install_vscode:

Visual Studio Code
------------------
Working with OpenFOAM involves a lot of editing of text files and running commands in a Linux shell. We do both in Visual Studio Code, which runs inside the container and which you open in your browser at http://localhost:8080. It opens the folder :code:`HydrothermalFoam_runs`.

* Open a terminal with *Terminal → New Terminal* (or CTRL + \`). This is a Linux shell inside the container, with OpenFOAM and HydrothermalFoam ready to use. All OpenFOAM commands in this class are typed in such a terminal.
* Open a file in the editor by clicking on it, or from the terminal with :code:`code filename`.
* Python files and Jupyter notebooks also work in Visual Studio Code. The course Python is already selected.

.. tip::
    You can also get a shell in the container from a terminal on your computer, without Visual Studio Code: :code:`docker compose exec htf bash`

JupyterLab
----------
JupyterLab runs inside the container as well; open it at http://localhost:8888. New notebooks use the course Python with all packages installed. Terminals in JupyterLab are Linux shells in the container, just like in Visual Studio Code.

Python
------
Python is already installed in the container with all packages we need in class: you do not have to install or activate anything. :code:`python` in a terminal, notebooks in JupyterLab, and Python files in Visual Studio Code all use the same course environment.

If you need an extra package, install it with :code:`pip install <package>` in a terminal in the container. Note that such packages are lost when the container is removed (see the warning above).

Test your installation
----------------------
Open a terminal in Visual Studio Code and check that OpenFOAM and HydrothermalFoam are available:

.. code-block:: bash

    echo $WM_PROJECT_VERSION    # should print 10
    which HTFoam

Then create a new notebook in JupyterLab (or a file :code:`hello.ipynb` in Visual Studio Code), and run this code:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import meshio

    x = np.linspace(0, 2 * np.pi, 100)
    plt.plot(x, np.sin(x))
    plt.show()

You should see a sine curve. If everything works, you are ready for the :ref:`L01_FirstCase`.

Troubleshooting
---------------

* **"Cannot connect to the Docker daemon"**: Docker Desktop is not running. Start it and try again.
* **"port is already allocated"**: another program uses port 8080 or 8888. Stop it, or change the first number in the :code:`ports` section of :code:`docker-compose.yml` (e.g. :code:`"127.0.0.1:8081:8080"`) and open http://localhost:8081 instead.
* **A new version of the course image**: if we announce an update, run this in your course folder (your files in :code:`HydrothermalFoam_runs` are kept):

  .. code-block:: bash

      docker compose pull
      docker compose up -d
