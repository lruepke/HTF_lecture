.. include:: /include.rst_

.. _L04_Exercise1:

Excercise 1
===========

We have prepared jupyter notebooks that let's you explore the thermodynamic properties of pure water. Play with them!

.. figure:: /_figures/Properties_Water.*
   :align: center
   :name: fig:Properties_Water

   Water properties as a function of temperature and pressure, using python package of iapws_. Properties calculation details can be found in the :ref:`/lectures/L04/cases/Jupp_Schultz/jupyter/Plot_Lectures.ipynb#Properties-as-function-of-temperature` in the notebook.

.. toctree::
    :maxdepth: 2

    cases/Jupp_Schultz/jupyter/Plot_Lectures.ipynb


.. admonition:: Assignment

   Try to get these notebooks to work and play a little bit with them. Explore how the values change under different P-T conditions and change the line plot to also show enthalpy. One way of doing this is to create a new notebook and to copy the code from above into them

   .. code-block:: bash 

      code plot_iapws_class.ipynb

.. tip:: Dangers and Annoyances

   Sometimes Visual Studio's way of setting the python interpretor is not that "intuitive". In principle, you set the the interpretor (your python environment) using CMD/CTRL+SHIFT+P and typing Python: Select Interpretor. Then you can chose the python environment of the class that has all the necessarry packages installed. 

   Problem is, .ipynb files also contain information on the interpretor that should be used with a notbook and Visual Studio uses that information. Therefore, if you download the scripts from us, there is a non-zero chance that your Visual Studio will chose the "wrong" interpretor and that you get an error message that certain packages are not installed. 

   In that case you need to follow a workaround. Rename your .ipynb file to .ipynb.json and open that file. Find the section on kernelspec:

   .. code-block:: bash

      "kernelspec": {
      "display_name": "Python 3",
      "language": "python",
      "name": "python3"
      }, 

   and delete that entire section. Save the file again as a .ipynb file and re-open it in Visual Studio - and hope for the best...
