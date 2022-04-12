Installation guide
==================

We will use Python for learning the basics of the finite element method. Later in the course, we will also use `FENICS <https://fenicsproject.org/>`_ to explore some more advanced concepts. Most of the work will be done in Jupyter notebooks. Let's get all of this to work.

Visual Studio Code
------------------
We will do a lot of editing of text files and you can use your favorite text editor for this. However, we  recommend to use `Microsoft's Visual Studio Code <https://code.visualstudio.com/>`_, which also nicely integrates with python. 

Python
--------
In case you already have a working python environment, you can adapt it for this course (e.g. by creating a new virtual environment). If not, we recommend `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Follow the miniconda installation instructions and afterwards create a virtual environment for this course. If you are asked to automatically activate the base environment (add it to the system path), chose "no". It's usually a good idea to keep the normal OS python environment intact and only activate a miniconda environment when you need it.

Download and install miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Just follow the installation instructions and keep the default options. We recommend to install python 3. During the installation, choose to **not** add it to your path (that's the default). Adding miniconda/anaconda to your :code:`$PATH` may seem convenient but there are several reason to not do it.

    * Python is used by many different tools on your computer, which probably expect that just calling python will use the Python (and additional packages) installed by the operating system. None of these will be available to Miniconda's Python.

    * The conda environment we will create contains several binary dependencies and we do not want to interfer with defaults on your system when, e.g. compiling software unrelated to our lecture.

Create a virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We will create a so-called virtual environment with all the python packages we will use during this class. To not infer with your default python installation, we will do this in a virtual environment. To get started open a terminal with activated base miniconda installation. 

.. admonition:: Starting python

    If you are on Windows, start an Anaconda Powershell Prompt from the start menu.

    On MacOS / Linux, open a terminal and type

    .. code-block:: bash

        conda activate base

    You can also do that in the terminal within Visual Studio Code (on MacOS).

Now we are ready to create a virtual environment. We can create it with this command:

.. code-block:: bash

    conda create -n py37_fem_class python=3.7 numpy pandas matplotlib vtk h5py ipython scipy ipykernel jupyterlab

We are using python 3.7 here (instead of the newest 3.8) because of an incompatibility with vtk. Activate the new environment

.. code-block:: bash

    conda activate py37_fem_class


Switching between environments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can activate and deactivate environments like this:

.. code-block:: bash

    conda activate py37_fem_class
    conda deactivate 

Working with jupyter notebooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will do most excercises using jupyter notebooks. A good workflow is to start jupyter labs in working directory. 

.. code-block:: bash

    cd "your working directory"
    jupyter lab


One possible issue is that you need to make sure that your jupyter notebooks use the correct python environment. One way of doing this is to register your virtual environment as a jupyter kernel. This can be done like this:

.. code-block:: bash

    conda activate py37_fem_class
    ipython kernel install --user --name=py37_fem_class

Restart you jupyter lab and try to select the correct python kernel.


Integration with Visual Studio Code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You will need to install Microsoft's Python extension. Just search for Python under Extensions and chose the one from Microsoft (usually the first option). Finally, you will have to set the Python interpreter. Do this by pushing CMD/CTRL+SHIFT+P. Type Python: Select Interpretor and select our newly created anaconda environment. If it doesn't show up, close and re-open Visual Studio Code.
