.. _Installation guide:

Installation guide
==================

We will use Python for learning the basics of the finite element method. Most of the work will be done in Jupyter notebooks. Let's get all of this to work.

Visual Studio Code
------------------
We will do a lot of editing of text files and you can use your favorite text editor for this. However, we  recommend to use `Microsoft's Visual Studio Code <https://code.visualstudio.com/>`_, which also nicely integrates with python. 

Python
--------

Download and install miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In case you already have a working python environment, you can adapt it for this course (e.g. by creating a new virtual environment). If not, we recommend `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Follow the miniconda installation instructions and afterwards create a virtual environment for this course. 

.. admonition:: conda-forge setup

    There are different sources, so-called channels, from where you can obtain the packages for your python environment. We will use the community channel *conda-forge*. 


Create a virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We will create a so-called virtual environment with all the python packages we will use during this class. To not interfere with your default python installation, we will do this in a virtual environment. To get started open a terminal with activated base miniconda installation. 

.. admonition:: Starting python

    If you are on Windows, start an Anaconda Powershell Prompt from the start menu.

    On MacOS / Linux, open a terminal and type

    .. code-block:: bash

        conda activate base

    You can also do that in the terminal within Visual Studio Code (on MacOS).


An easy way to create the environment is to use an :download:`environment.yaml <_static/files/environment.yaml>` file. This file contains all the necessary information to create the environment. Just copy it into a file and save it as environment.yaml.

.. code-block:: yaml

    name: py312_fem_class

    channels:
    - conda-forge
    - nodefaults

    dependencies:
    - python=3.12
    - ipykernel
    - matplotlib
    - numpy
    - pandas
    - vtk
    - pyevtk
    - h5py
    - scipy
    - ipykernel
    - pip

    - pip:
        - meshio
        - triangle


.. admonition:: Python versions

    In Summer 2025, conda already supports python 3.13. However, some packages are not working with the newest python and we therefore use python 3.12. Please keep this in mind, if you are working with your own python installation.


Now create the environment by typing

.. code-block:: bash

    conda env create -f environment.yaml

If everything works, you should see something like this at the end of the output:

.. code-block:: bash

    #
    # To activate this environment, use
    #
    #     conda activate py312_fem_class
    #
    # To deactivate an active environment, use
    #
    #     conda deactivate
    #

Activate the environment
^^^^^^^^^^^^^^^^^^^^^^^^^^
You can now activate the environment by typing

.. code-block:: bash

    conda activate py312_fem_class


Switching between environments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can activate and deactivate environments like this:

.. code-block:: bash

    conda activate py312_fem_class
    conda deactivate 


Integration with Visual Studio Code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You will need to install Microsoft's Python extension. Just search for Python under Extensions and chose the one from Microsoft (usually the first option). Finally, you will have to set the Python interpreter. Do this by pushing CMD/CTRL+SHIFT+P. Type Python: Select Interpreter and select our newly created anaconda environment. If it doesn't show up, close and re-open Visual Studio Code.

We prefer running the jupyter notebooks in Visual Studio Code but you can also use jupyter lab.
