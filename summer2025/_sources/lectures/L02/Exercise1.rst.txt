.. include:: /include.rst_

.. _L02_Exercise1:


Excercise 1
=======================================

Transient flow in a lid-driven cavity
--------------------------------------

In this section, we explore how the fluid flow in a lid-driven cavity develops over time, and how the *behavoir* and *time scales* of the simulation are influenced by model parameters. Understanding this helps to interpret simulation results and perform meaningful parametric studies.

Time scales and regimes of viscous flow
------------------------------------------

Reynolds number
^^^^^^^^^^^^^^^

One important dimensionless number in fluid dynamics is the **Reynolds number**. It is defined as:

.. math::
    :label: eq:Re_cavity
    
    Re = \frac{UL}{\nu}


where:

- :math:`U` is the characteristic velocity (e.g., lid speed),
- :math:`L` is the characteristic length scale of the domain (e.g., cavity height or width),
- :math:`\nu` is the kinematic viscosity.


The Reynolds number indicates the relative importance of inertial forces (which cause the fluid to keep moving) to viscous forces (which resist motion) in the flow. In the case of a lid-driven cavity, a low Reynolds number (e.g., :math:`Re < 1`) indicates that viscous forces dominate, leading to laminar flow. A moderate Reynolds number of  :math:`100 < Re < 2000` indicates that inertial forces start to compete with viscous forces, and the flow may still be laminar but with more complex patterns. A high Reynolds number (e.g., :math:`Re > 2000`) indicates that inertial forces dominate, leading to turbulent flow.


Let's carefully check our base case. The spatial dimensions are set in :code:`system/blockMeshDict`:


.. code-block:: foam 
    :name: lst:2dcavity_blockMeshDict
    :linenos:
    :emphasize-lines: 12, 20, 21, 33
    :caption: blockMeshDict of the cavity flow tutorial. 

    /*--------------------------------*- C++ -*----------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     | Website:  https://openfoam.org
      \\  /    A nd           | Version:  9
       \\/     M anipulation  |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        format      ascii;
        class       dictionary;
        object      blockMeshDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    convertToMeters 0.1;

    vertices
    (
        (0 0 0)
        (1 0 0)
        (1 1 0)
        (0 1 0)
        (0 0 0.1)
        (1 0 0.1)
        (1 1 0.1)
        (0 1 0.1)
    );


    blocks
    (
        hex (0 1 2 3 4 5 6 7) (20 20 1) simpleGrading (1 1 1)
    );



We are defining a 3-D mesh with a width and height of 0.1 m and a thickness of 0.01 m.  A potential pitfall is line 12, where a global scaling of the spatial dimensions is introduced. The 3rd dimension does not matter here, as we are running the simulation in 2d mode (boundary faces are set to ``empty``). The lid velocity is set in :code:`0/U` to 1 m/s. The kinematic viscosity is set in :code:`constant/transportProperties` to 0.01 m^2/s. The numerical resolution is set in :code:`system/blockMeshDict` to 20 cells in the x and y direction, which gives us a grid spacing of 0.005 m. 

The Reynolds number for this case is:

.. math::
    :label: eq:Re_cavity_num

    Re = \frac{UL}{\nu} = \frac{1~\text{m/s} \cdot 0.1~\text{m}}{0.01~\text{m}^2/\text{s}} = 10


This indicates that the flow is laminar, as the Reynolds number is below 1000.


Time scales
^^^^^^^^^^^^

 Let's check the time scales of the problem. One time scale of interest is the **convection time scale**, which is given by:

.. math::
    :label: eq:t_conv

    t_\text{conv} \sim \frac{L}{U}


where:

- :math:`L` is the characteristic length scale of the domain (e.g., cavity height),
- :math:`U` is the characteristic velocity (e.g., lid speed).


This relation tells us how long it takes for a fluid parcel to travel the length of the domain. For our case, this is:

.. math::
    :label: eq:t_conv_num

    t_\text{conv} \sim \frac{0.1~\text{m}}{1~\text{m/s}} = 0.1~\text{s}


Another important time scale is the **viscous diffusion time scale**, which is given by:

.. math::
    :label: eq:t_visc

    t_\text{visc} \sim \frac{L^2}{\nu}


This relation tells us how quickly momentum spreads through the fluid. Higher viscosity or smaller domains lead to faster development of flow patterns. For our case, this is:

.. math::
    :label: eq:t_visc_num

    t_\text{visc} \sim \frac{(0.1~\text{m})^2}{0.01~\text{m}^2/\text{s}} = 1~\text{s}


With these time scales, we can see that the convection time scale is  shorter than the viscous diffusion time scale. This means that the flow will develop quickly and reach a steady state in a short time. That's why we don't see any significant changes in the flow pattern in our base case, which runs at :math:`\delta t = 0.005 s` until 0.5 s, outputting results every 0.1 s.


Courant number and numerical stability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, we need to check the **Courant number**, which is a measure of the stability of the numerical scheme (we will get to the details later). The Courant number is defined as:

.. math::
    :label: eq:Co

    Co = \frac{U \Delta t}{\Delta x}


where:

- :math:`U` is the characteristic velocity (e.g., lid speed),
- :math:`\Delta t` is the time step size,
- :math:`\Delta x` is the grid spacing (e.g., cell size).


The Courant number indicates how far a fluid parcel moves in one time step relative to the grid spacing. A Courant number less than 1 indicates that the numerical scheme is stable, while a Courant number greater than 1 indicates that the scheme may become unstable.
For our case, the Courant number is:

.. math::
    :label: eq:Co_num

    Co = \frac{1~\text{m/s} \cdot 0.005~\text{s}}{0.005~\text{m}} = 1


This indicates that the numerical scheme is stable (just about), as the Courant number is less than 1. However, if we were to increase the lid speed or decrease the grid spacing, we would need to decrease the time step size to maintain a stable Courant number.



Exploration Framework
---------------------

Let's modify the base case and explore how the flow develops under different physical and numerical conditions. 


1. **Baseline Case**
   - Use the standard lid-driven cavity example with:
     - Domain size: :math:`L = 0.1~\text{m}`
     - Viscosity: :math:`\nu = 0.01~\text{m}^2/\text{s}`
     - Lid speed: :math:`U = 1~\text{m/s}`
   - Observe flow development over time.
   - Questions to explore:
     - How long does it take for the flow to reach a steady state? Change the output interval and total run time to check the transient flow development.
     - What is the shape of the vortex?
     - How does the velocity field look like?

2. **Vary Viscosity**
   - Change :math:`\nu` to :math:`1 \times 10^{-1}` and :math:`1 \times 10^{-3}`.
   - Questions to explore:
     - How fast does the vortex form?
     - How does the flow pattern change?

3. **Change Lid Velocity**
   - Keep :math:`\nu` fixed and vary lid velocity :math:`U`.
   - Questions to explore:
     - Does a higher velocity lead to a longer or shorter transient?
     - Does the vortex strength or location change?

5. **Mesh Sensitivity**
   - Increase mesh resolution (e.g., from 20x20 to 40x40).
   - Observe differences in vortex shape and location.
   - Discuss numerical vs physical convergence.


Expected Observations
---------------------

- **High viscosity**: Flow develops quickly and becomes smooth and laminar.
- **Low viscosity**: Transient lasts longer, possible oscillations or instabilities.
- **High lid velocity**: Stronger vortices, higher Reynolds number, more complex flow.
- **Fine mesh**: Captures finer details (e.g., secondary vortices), more accurate but slower to compute.

