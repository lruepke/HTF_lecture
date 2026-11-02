.. include:: /include.rst_

.. _L02_Exercise1:


Exercise 1
=======================================

Transient flow in a lid-driven cavity
--------------------------------------

In this section, we explore how the fluid flow in a lid-driven cavity develops over time, and how the *behavior* and *time scales* of the simulation are influenced by model parameters. Understanding this helps to interpret simulation results and perform meaningful parametric studies.

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


The Reynolds number compares inertia (moving fluid tends to keep moving) with viscosity (friction within the fluid, which slows the motion down):

- **Small Re**: viscosity wins. The flow is smooth and predictable, and the fluid moves in orderly layers. This is called *laminar* flow.
- **Large Re**: inertia wins. The flow becomes more complex, with additional vortices and eventually unsteady flow, and at very large Re it becomes *turbulent* (chaotic).

The values of Re at which these changes happen depend on the geometry of the problem, so there is no single threshold that applies to all flows.

.. admonition:: Why this matters for hydrothermal flow

    Later in the course we model fluid flow through porous rocks. There, the pores are tiny and the flow is very slow, so the Reynolds number is much smaller than 1. Viscosity wins completely, and the flow can be described by the much simpler Darcy's law instead of the full Navier-Stokes equations.


Let's carefully check our base case. The spatial dimensions are set in :code:`system/blockMeshDict`:


.. code-block:: foam 
    :name: lst:2dcavity_blockMeshDict
    :linenos:
    :emphasize-lines: 16, 21, 22, 33
    :caption: blockMeshDict of the cavity flow tutorial. 

    /*--------------------------------*- C++ -*----------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     | Website:  https://openfoam.org
      \\  /    A nd           | Version:  10
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



We are defining a 3-D mesh with a width and height of 0.1 m and a thickness of 0.01 m.  A potential pitfall is line 16, where a global scaling of the spatial dimensions is introduced. The 3rd dimension does not matter here, as we are running the simulation in 2d mode (boundary faces are set to ``empty``). The lid velocity is set in :code:`0/U` to 1 m/s. The kinematic viscosity is set in :code:`constant/physicalProperties` to 0.01 m^2/s. The numerical resolution is set in :code:`system/blockMeshDict` to 20 cells in the x and y direction, which gives us a grid spacing of 0.005 m. 

The Reynolds number for this case is:

.. math::
    :label: eq:Re_cavity_num

    Re = \frac{UL}{\nu} = \frac{1~\text{m/s} \cdot 0.1~\text{m}}{0.01~\text{m}^2/\text{s}} = 10


This is a small Reynolds number: viscosity dominates, and we expect a smooth, laminar flow with a single large vortex.

To get a feeling for the numbers, here are some flows from the Earth sciences. The kinematic viscosity is the dynamic viscosity :math:`\eta` (in Pa s) divided by the density: :math:`\nu = \eta/\rho`.

.. list-table:: Reynolds numbers of some geological flows (rough values).
   :header-rows: 1
   :widths: 28 18 16 20 18

   * - Flow
     - Velocity :math:`U`
     - Length :math:`L`
     - Viscosity :math:`\nu`
     - :math:`Re`
   * - Mantle flow below a moving tectonic plate
     - 5 cm/yr :math:`\approx 2 \times 10^{-9}` m/s
     - 1000 km
     - :math:`10^{21}` Pa s / 3300 kg/m³ :math:`\approx 3 \times 10^{17}` m²/s
     - :math:`\approx 10^{-20}`
   * - Hot water flowing through the pores of a rock
     - 1 m/yr :math:`\approx 3 \times 10^{-8}` m/s
     - 1 mm (pore size)
     - :math:`10^{-7}` m²/s
     - :math:`\approx 10^{-6}`
   * - Our cavity
     - 1 m/s
     - 0.1 m
     - 0.01 m²/s
     - 10
   * - An ocean current like the Gulf Stream
     - 1 m/s
     - 100 km
     - :math:`10^{-6}` m²/s (water)
     - :math:`\approx 10^{11}`

The range is enormous. The mantle is so viscous that inertia plays no role at all: if a plate stopped moving, the mantle flow below it would stop practically at once. This is why mantle convection is modeled with the Stokes equations, which are the Navier-Stokes equations without the inertia terms. The ocean is the other extreme: its currents are fully turbulent and full of eddies of all sizes. The flow through the pores of a rock, which is what this course is about, is again dominated by viscosity (see the box "Why this matters for hydrothermal flow" above).


Time scales
^^^^^^^^^^^^

 Let's check the time scales of the problem. One time scale of interest is the **convection time scale**, which is given by:

.. math::
    :label: eq:t_conv

    t_\text{conv} \sim \frac{L}{U}


where:

- :math:`L` is the characteristic length scale of the domain (e.g., cavity height or length),
- :math:`U` is the characteristic velocity (e.g., lid speed).


This relation tells us how long it takes for a fluid parcel to travel the length of the domain. For our case, this is:

.. math::
    :label: eq:t_conv_num

    t_\text{conv} \sim \frac{0.1~\text{m}}{1~\text{m/s}} = 0.1~\text{s}


Another important time scale is the **viscous diffusion time scale**, which is given by:

.. math::
    :label: eq:t_visc

    t_\text{visc} \sim \frac{L^2}{\nu}


This relation tells us how quickly the motion of the lid spreads into the fluid. Higher viscosity or smaller domains lead to faster development of flow patterns. For our case, this is:

.. math::
    :label: eq:t_visc_num

    t_\text{visc} \sim \frac{(0.1~\text{m})^2}{0.01~\text{m}^2/\text{s}} = 1~\text{s}

.. admonition:: You know this from heat conduction

    The time it takes heat to diffuse over a distance :math:`L` is :math:`t \sim L^2/\kappa`, with the thermal diffusivity :math:`\kappa`. This is how we estimate how long a dike or the oceanic lithosphere needs to cool. Viscosity works the same way: the kinematic viscosity :math:`\nu` is the "diffusivity of momentum", and it has the same units as :math:`\kappa` (m²/s). The fluid in the cavity starts at rest, and the motion of the lid diffuses into it just like heat diffuses into cold rock.

So :math:`t_\text{visc}` is our estimate for how long it takes until the whole cavity is moving and the flow is steady. Such estimates tell you whether something takes milliseconds, seconds or hours. They are not exact: expect them to be off by a factor of ten or more. How good is the estimate here? Our base case runs until 0.5 s and writes results every 0.1 s. You will measure the time to steady state yourself in the first task below.


Courant number and numerical stability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, we need to check the **Courant number**. It tells us whether our time step is small enough for the mesh we use (we will get to the details later). The Courant number is defined as:

.. math::
    :label: eq:Co

    Co = \frac{U \Delta t}{\Delta x}


where:

- :math:`U` is the characteristic velocity (e.g., lid speed),
- :math:`\Delta t` is the time step size,
- :math:`\Delta x` is the grid spacing (e.g., cell size).


The Courant number tells us how many cells the fluid moves in one time step. The rule of thumb is: **in one time step, the fluid should not move further than one cell**, so :math:`Co` should stay below 1. If it gets larger, the simulation can become unstable: the results start to oscillate and the run eventually crashes.
For our case, the Courant number is:

.. math::
    :label: eq:Co_num

    Co = \frac{1~\text{m/s} \cdot 0.005~\text{s}}{0.005~\text{m}} = 1


This is just about acceptable. If we increase the lid speed or make the cells smaller, we have to make the time step smaller as well.

You do not have to calculate the Courant number by hand: the solver reports it in every time step. If you run the case with :code:`./run.sh`, look into the log file:

.. code-block:: bash

    grep "Courant" log.icoFoam | tail -1

.. code-block:: none

    Courant Number mean: 0.222158 max: 0.852134

The maximum is a bit smaller than our estimate, because the fluid inside the cavity moves slower than the lid. Check this line whenever you change the velocity, the mesh, or the time step.



Exploration Framework
---------------------

Let's modify the base case and explore how the flow develops under different physical and numerical conditions. 


1. **Baseline Case: predict, then measure**

   - Use the standard lid-driven cavity example with:

     - Domain size: :math:`L = 0.1~\text{m}`
     - Viscosity: :math:`\nu = 0.01~\text{m}^2/\text{s}`
     - Lid speed: :math:`U = 1~\text{m/s}`

   - **Predict first:** before you run anything, write down when you expect the flow to be steady, based on the time scales above.

   - **Measure:** set :code:`writeInterval 2;` in :code:`system/controlDict`, so that results are written every 0.01 s, and run the case again. In ParaView, select :code:`a.foam` and apply *Filters → Data Analysis → Probe Location* with the point in the center of the cavity, (0.05 0.05 0.005). Then apply *Filters → Data Analysis → Plot Data Over Time* to the probe and plot the velocity magnitude against time.

   - Questions to explore:

     - After how much time does the velocity in the center stop changing? Compare this with your prediction. By what factor were you off, and why can a scaling estimate be off by that much?
     - Look at the time steps written by the original case (every 0.1 s). Why did the original case look so "boring"?
     - What is the shape of the vortex, and where is its center?
     - How does the velocity field look like?

2. **Vary Viscosity**
   
   - Change the viscosity in :code:`constant/physicalProperties` from :math:`\nu = 0.01~\text{m}^2/\text{s}` to :math:`0.1~\text{m}^2/\text{s}` and to :math:`0.001~\text{m}^2/\text{s}`. Use a copy of the case for each value.

   - Before you run: calculate :math:`Re` and :math:`t_\text{visc}` for each value. Does the run time of 0.5 s still make sense? For :math:`\nu = 0.001~\text{m}^2/\text{s}`, increase :code:`endTime` in :code:`system/controlDict` (e.g. to 2 s).

   - Questions to explore:

     - How fast does the vortex form? Measure it as in task 1 and compare with your estimate. Is the trend right?
     - How does the flow pattern change? Where is the center of the vortex?

3. **Change Lid Velocity**
   
   - Keep :math:`\nu` fixed and vary lid velocity :math:`U`.

   - Questions to explore:

     - Does a higher velocity lead to a longer or shorter transient?
     - Does the vortex strength or location change?

4. **Mesh Sensitivity**
   
   - Increase mesh resolution (e.g., from 20x20 to 40x40).

   - Observe differences in vortex shape and location.

   - Discuss numerical vs physical convergence.



.. admonition:: Things to keep in mind

        - You may have to adjust the time step size to maintain a stable Courant number. 
        - When changing time step size, make sure to adjust the output interval accordingly.
        - Make copies of the case for each variation to keep track of different simulations.
        - You can visualize multiple simulations in ParaView by loading multiple case folders at once.
        - Use the "plot over line" filter in ParaView to extract velocity profiles and compare different cases
        - Use the "stream tracer" filter in ParaView to visualize the flow field and vortex structure.




Expected Observations
---------------------

- **Baseline**: the velocity in the center reaches 99 % of its final value after about 0.07 s, more than ten times faster than the estimate :math:`t_\text{visc} \sim 1` s. So the original case was already steady at its first output time (0.1 s). The estimate gets the trend right (twice the viscosity, twice as fast), but not the exact number. See the box below for why.
- **High viscosity**: Flow develops quickly and becomes smooth and laminar.
- **Low viscosity**: Transient lasts longer, possible oscillations or instabilities.
- **High lid velocity**: Stronger vortices, higher Reynolds number, more complex flow.
- **Fine mesh**: Captures finer details (e.g., secondary vortices), more accurate but slower to compute.

.. admonition:: Why is the flow steady after 0.07 s and not after 1 s?

    The estimate :math:`t_\text{visc} \sim L^2/\nu` only tells us the order of magnitude. Two things make the real time shorter:

    1. The motion of the lid does not have to diffuse across the whole cavity. The flow pattern is set up once the motion has reached about the center of the cavity, which is only half the distance, :math:`L/2`. Because the distance enters squared, half the distance needs only a quarter of the time.
    2. The distance over which something diffuses in the time :math:`t` is about :math:`2\sqrt{\nu t}`, not :math:`\sqrt{\nu t}`. You know this from the cooling of the oceanic lithosphere, whose thickness grows as about :math:`2.3\sqrt{\kappa t}`.

    Together, the motion reaches the center when :math:`2\sqrt{\nu t} \approx L/2`, i.e. after

    .. math::

        t \approx \frac{L^2}{16\,\nu} = \frac{(0.1~\text{m})^2}{16 \cdot 0.01~\text{m}^2/\text{s}} \approx 0.06~\text{s},

    close to what we measured. The dependence on :math:`\nu` stays the same: with a ten times smaller viscosity, :math:`\nu = 0.001~\text{m}^2/\text{s}`, the flow needs about 0.7 s, ten times longer. The lesson: scaling estimates are great for trends and for comparing cases, but numerical factors like :math:`1/16` can make a big difference in absolute numbers.

