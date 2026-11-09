.. include:: /include.rst_

.. _L02_Debugging:

Exercise: Fix a broken case
===========================

Sooner or later, every OpenFOAM model stops with an error message. This is normal, even for experienced users. The error messages look scary at first, but they usually tell you quite precisely what is wrong and where. In this exercise, you practice reading them.

Getting the cases
-----------------

Download :download:`broken_cases.zip <cases/broken_cases.zip>`, move it into the folder :code:`HydrothermalFoam_runs` in your course folder, and unpack it in a terminal in the container:

.. code-block:: bash

    cd $HOME/HydrothermalFoam_runs
    python -m zipfile -e broken_cases.zip .
    cd broken_cases
    ls

There are four cases, :code:`broken1` to :code:`broken4`. Each is the cavity case from this lecture with **one** error in it.

Your task
---------

For each case:

1. Go into the case directory and run the mesh generator and the solver directly, so that you see the error message in the terminal:

   .. code-block:: bash

       cd broken1
       blockMesh
       icoFoam

2. Find the block that starts with :code:`--> FOAM FATAL ERROR` or :code:`--> FOAM FATAL IO ERROR` and read it carefully. Does it name a case file, a keyword, or a patch?
3. Open the file in Visual Studio Code, find the error and fix it.
4. Run :code:`icoFoam` again (and :code:`blockMesh`, if you changed the mesh) until the case runs to the end. Then go to the next case with :code:`cd ../broken2`.

For each case, write down what the error message said and what the actual error was. Did the message point directly to the problem?

.. admonition:: How to read an OpenFOAM error message

    .. code-block:: none

        --> FOAM FATAL IO ERROR:
        Cannot find patchField entry for movingWall

        file: .../0/p/boundaryField from line 24 to line 34.

            From function void Foam::GeometricBoundaryField<...>::readField(...)
            in file .../GeometricBoundaryField.C at line 187.

        FOAM exiting

    * The first lines say **what** is wrong.
    * The line starting with :code:`file:` says **where**: the file in your case, and often the line numbers.
    * The lines starting with :code:`From function` and :code:`in file` refer to the OpenFOAM source code. You can usually ignore them.

.. tip::
    * The error is not always exactly where the message points. If OpenFOAM complains about a keyword that is clearly there, look at the line **before** it.
    * A case that runs without an error message is not necessarily correct. Always check the log, e.g. the Courant number, and look at the results.
    * If you are stuck, compare your case with the original tutorial. :code:`diff` shows all lines that differ:

      .. code-block:: bash

          diff -r $FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity .

      Do not peek too early, though: reading the error messages is the point of this exercise!
