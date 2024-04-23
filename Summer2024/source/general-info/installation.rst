Installation guide
==================

We will use Python 3 for learning the basics of the finite element method. Later in the course, we will also use `FENICS <https://fenicsproject.org/>`_ to explore some more advanced concepts. Most of the work will be done in Jupyter notebooks. Let's get all of this to work.

Visual Studio Code
------------------
We will do a lot of editing of text files and you can use your favorite text editor for this. However, we  recommend to use `Microsoft's Visual Studio Code <https://code.visualstudio.com/>`_, which also nicely integrates with python. 

Python
--------

Download and install miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In case you already have a working python environment, you can adapt it for this course (e.g. by creating a new virtual environment). If not, we recommend `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Follow the miniconda installation instructions and afterwards create a virtual environment for this course. If you are asked to automatically activate the base environment (add it to the system path), chose "no". It's usually a good idea to keep the normal OS python environment intact and only activate a miniconda environment when you need it.

.. admonition:: conda-forge setup

    There are different sources, so-called channels, from where you can obtain the packages for your python environment. We will use the community channel *conda-forge*. 

    Our recommended setup is to the use a very basic *base environment* that only cotains the necessarry packages to run *jupyter notebook* and/or *jupyter lab*. Everything else will be done from virtual enviroments, which you can activate from within jupyter. 

    .. code-block:: bash

        conda install -n base -c conda-forge jupyterlab notebook nb_conda_kernels


Create a virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We will create a virtual environment with all the python packages we will use during this class. To not interfer with your default python installation, we will do this in a virtual environment. To get started open a terminal with activated base miniconda installation. 

Make sure cour conda base environment is activated

.. code-block:: bash

    conda activate base

Now we are ready to create a virtual environment. We can create it with this command:

.. code-block:: bash

    conda create -n py3_fem_class numpy pandas matplotlib scipy ipykernel

Activate the new environment

.. code-block:: bash

    conda activate py3_fem_class


Switching between environments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can activate and deactivate environments like this:

.. code-block:: bash

    conda activate py3_fem_class
    conda deactivate 

Working with jupyter notebooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will do most excercises using jupyter notebooks. A good workflow is to start jupyter notebook in the conda base environment and then chose the the python kernel (the virtual environment) inside the notebook.

.. code-block:: bash

    cd "your working directory"
    conda activate base
    jupyter notebook


Now create a new notebook and choose *py3_fem_class* as your kernel. Check that you can import e.g. pandas. 
 

.. admonition:: Confused?

    If you have never used python or are new to jupyter notebooks, no worries! Things will become clear when we are doing the actual exercises. 

Integration with Visual Studio Code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You will need to install Microsoft's Python extension. Just search for Python under Extensions and chose the one from Microsoft (usually the first option). Finally, you will have to set the Python interpreter. Do this by pushing CMD/CTRL+SHIFT+P. Type Python: Select Interpretor and select our newly created anaconda environment. If it doesn't show up, close and re-open Visual Studio Code.
