.. include:: /include.rst_

.. _L05_Exercise1:

Exercise 1: Fluid properties
============================

Let's start by computing some properties of pure water as a function of pressure and temperature.

.. figure:: /_figures/density_pt.*
   :align: center
   :name: fig:Properties_Water

   Pure water density as a function of temperature and pressure. Plot was made using the python package of iapws_.



We can use the python package iapws to compute the fluid properties. One way to do this is to call the IAPWS97 class and to specify the temperature and pressure. The class then returns a build-in datatype that holds all the properties.

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   import iapws


.. code-block:: python

   # pressure and temperature in MPa, and C
   p, T =35, 300

   #compute properties
   steam=iapws.IAPWS97(T=T+273.15, P=p) # note that iaapws uses temperature in K and pressure in MPa

   #format the output for printing
   output_format = "{:<17} = {:>10.2f} {}\n"
   viscosity_format = "{:<17} = {:>10.2e} {}\n"

   output = (output_format.format("Density", steam.rho, "kg/m\u00B3") +
            output_format.format("Specific heat", steam.cp*1000, "J/(kg\u00B7K)") +
            viscosity_format.format("Viscosity", steam.mu, "Pa\u00B7s") +
            output_format.format("Specific enthalpy", steam.h, "J/kg"))

   print(output)# import water properties

.. code-block:: bash

   Density           =     757.72 kg/m³
   Specific heat     =    4990.59 J/(kg·K)
   Viscosity         =   9.47e-05 Pa·s
   Specific enthalpy =    1326.81 J/kg


Fluid properties as function of p and T
------------------------------------------

Next we want to make plots of the fluid properties as a function of pressure and temperature. We can do this by looping over the pressure and temperature and computing the properties for each combination. We then store the properties and plot them.

Take this notebook as a starting point and make nice looking plots for the other fluid properties such as viscosity, specific heat, and enthalpy. You can find the documentation of the iapws package here: https://iapws.readthedocs.io/en/latest/

.. toctree::
    :maxdepth: 2

    plot_water_properties.ipynb


What happens when you chose a pressure below the critical pressure at 22.064 MPa?

Fluxibility
------------------------------------------

Use the same notbook to make plots of fluxibility as a function of temperature for different pressures. Would we expect different vent temperatures as a function of depth of the driving heat source?


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

   and delete that entire section. Save the file again as a .ipynb file and re-open it in Visual Studio - and hope for the best.
