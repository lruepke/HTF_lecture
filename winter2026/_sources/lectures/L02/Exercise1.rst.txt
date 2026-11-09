.. include:: /include.rst_

.. _L02_Exercise1:


Exercise 1
=======================================

Transient flow in a lid-driven cavity
--------------------------------------

In this section, we explore how the fluid flow in a lid-driven cavity develops over time, and how the *flow pattern* and *time scales* are influenced by model parameters. Understanding this helps to interpret simulation results and perform meaningful parametric studies.

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

- **Small Re**: viscous forces dominate inertia. We expect smooth, *laminar* flow, in which the fluid moves in orderly layers.
- **Large Re**: inertia becomes more important. Depending on the geometry and forcing, the flow may develop additional vortices, become unsteady, and eventually become *turbulent*, i.e. chaotic.

The values of Re at which these changes happen depend on the geometry of the problem, so there is no single threshold that applies to all flows.

.. admonition:: Why this matters for hydrothermal flow

    Later in the course we model fluid flow through porous rocks. For the porous-flow cases we study, pore sizes and flow speeds are small, so the Reynolds number based on pore size and pore-scale velocity is much smaller than 1. Viscosity wins completely, and the flow can be described by the much simpler Darcy's law instead of the full Navier-Stokes equations.


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



We are defining a 3-D mesh with a width and height of 0.1 m and a thickness of 0.01 m.  A potential pitfall is line 16, where a global scaling of the spatial dimensions is introduced. The mesh is one cell thick in the z direction, with ``empty`` conditions on the front and back faces, so the solution is two-dimensional. The lid velocity is set in :code:`0/U` to 1 m/s. The kinematic viscosity is set in :code:`constant/physicalProperties` to :math:`0.01~\text{m}^2/\text{s}`. The numerical resolution is set in :code:`system/blockMeshDict` to 20 cells in each of the x and y directions, which gives us a grid spacing of 0.005 m.

The Reynolds number for this case is:

.. math::
    :label: eq:Re_cavity_num

    Re = \frac{UL}{\nu} = \frac{1~\text{m/s} \cdot 0.1~\text{m}}{0.01~\text{m}^2/\text{s}} = 10


This is a small Reynolds number: viscosity dominates, and we expect a smooth, laminar flow with one dominant vortex. Much weaker corner vortices may also be present.

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

The range is enormous. The mantle is so viscous that inertia plays no role at all: if the driving forces vanished, inertia would not keep the mantle moving for an appreciable time. This is why mantle convection is modeled with the Stokes equations, which are the Navier-Stokes equations without the inertia terms. The ocean is the other extreme: its currents are fully turbulent and full of eddies of all sizes. The flow through the pores of a rock, which is what this course is about, is again dominated by viscosity (see the box "Why this matters for hydrothermal flow" above).


Time scales
^^^^^^^^^^^^

Let's check the time scales of the problem. One time scale of interest is the **convection time scale**, which is given by:

.. math::
    :label: eq:t_conv

    t_\text{conv} \sim \frac{L}{U}


where:

- :math:`L` is the characteristic length scale of the domain (e.g., cavity height or length),
- :math:`U` is the characteristic velocity (e.g., lid speed).


This estimates the time to travel a distance :math:`L` at speed :math:`U`. Fluid inside the cavity moves more slowly than the lid, so this is not the actual circulation time of a fluid parcel. For our case, this is:

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

So :math:`t_\text{visc}` gives a characteristic time for momentum to diffuse across the cavity. It is a useful starting estimate for the transient duration, rather than a prediction of when every part of the flow reaches steady state. Such estimates tell you whether something takes milliseconds, seconds or hours. They are not exact: expect them to be off by a factor of ten or more. How good is the estimate here? Our base case runs until 0.5 s and writes results every 0.1 s. You will measure the time to steady state yourself in the first task below.


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


The Courant number tells us how many cells the fluid moves in one time step. The rule of thumb is: **in one time step, the fluid should not move further than one cell**, so :math:`Co` should stay below 1. If it gets larger, the results become less accurate, and the simulation can become unstable: the values start to oscillate and grow until they are meaningless.
For our case, the Courant number is:

.. math::
    :label: eq:Co_num

    Co = \frac{1~\text{m/s} \cdot 0.005~\text{s}}{0.005~\text{m}} = 1


This estimate is at the upper limit of our guideline. Check the maximum Courant number reported by the solver and reduce the time step if needed. If we increase the lid speed or make the cells smaller, we have to make the time step smaller as well.

You do not have to calculate the Courant number by hand: the solver reports it in every time step. If you run the case with :code:`./run.sh`, look into the log file:

.. code-block:: bash

    grep "Courant" log.icoFoam | tail -1

.. code-block:: none

    Courant Number mean: 0.222158 max: 0.852134

The reported maximum can be smaller than the estimate based on lid speed: OpenFOAM computes it from face fluxes, and the interior fluid generally moves more slowly than the lid. Check the maximum throughout the run whenever you change the velocity, mesh, or time step. A Courant number below 1 is a practical guideline for :code:`icoFoam`, not a guarantee of accuracy; time-step sensitivity also matters.



Exploration Framework
---------------------

Let's modify the base case and explore how the flow develops under different physical and numerical conditions. Most of the tasks below compare several runs, e.g. with different viscosities. Here is how to set this up.

.. admonition:: Comparing runs side by side

    **One folder per run.** Never change a case that you want to compare with. Make a copy for every variation and give it a name that says what is different, e.g. for a smaller viscosity:

    .. code-block:: bash

        cd $HOME/HydrothermalFoam_runs
        cp -r cavity2D cavity2D_nu0.001
        cd cavity2D_nu0.001

    Change the parameter, e.g. :code:`nu` in :code:`constant/physicalProperties`, and run the case with :code:`./run.sh`. The script cleans the copied results first.

    **Results next to each other in ParaView.**

    1. Open the :code:`a.foam` files of both cases (*File → Open*, once for each case folder) and click *Apply*. Both cases now lie on top of each other.
    2. Select the second case in the *Pipeline Browser* and apply *Filters → Alphabetical → Transform*. In its properties, set *Translate* to :code:`0.12 0 0` and click *Apply*. This moves the second case 0.12 m to the right, next to the first one. For a third case, use :code:`0.24 0 0`, and so on.
    3. If the original second case is still visible, hide it with the eye icon in the *Pipeline Browser*, so that only the moved copy is shown.
    4. Color both by the same field and give them the **same color range** (*Rescale to Custom Range*). Otherwise, the same color means different values in the two cases.

    Further filters, e.g. a *Calculator* or *Surface LIC*, are applied to the transformed case. To compare profiles, apply *Plot Over Line* to the **original** cases, not to the moved one, so that you can use the same line for both. Show both lines in the same *Line Chart View*: then they are directly on top of each other.


1. **Baseline Case: predict, then measure**

   - Use the standard lid-driven cavity example with:

     - Domain size: :math:`L = 0.1~\text{m}`
     - Viscosity: :math:`\nu = 0.01~\text{m}^2/\text{s}`
     - Lid speed: :math:`U = 1~\text{m/s}`

   - **Predict first:** before you run anything, write down when you expect the flow to be steady, based on the time scales above.

   - **Measure:** set :code:`writeInterval 2;` in :code:`system/controlDict`, so that results are written every 0.01 s, and run the case again. In ParaView, select :code:`a.foam` and apply *Filters → Data Analysis → Probe Location* with the point in the center of the cavity, (0.05 0.05 0.005). Then apply *Filters → Data Analysis → Plot Data Over Time* to the probe and plot the velocity magnitude against time. Use the final saved value as an initial reference: when does the magnitude stay within 1 % of it? Extend the run if the final value is still changing. This measures settling at one point; check another point or a velocity profile before calling the whole flow steady.

   - Questions to explore:

     - After how much time does the velocity in the center stop changing? Compare this with your prediction. By what factor were you off, and why can a scaling estimate be off by that much?
     - Look at the time steps written by the original case (every 0.1 s). Why did the original case look so "boring"?
     - What is the shape of the vortex, and where is its center?
     - What does the velocity field look like?

2. **Vary Viscosity**
   
   - Change the viscosity in :code:`constant/physicalProperties` from :math:`\nu = 0.01~\text{m}^2/\text{s}` to :math:`0.1~\text{m}^2/\text{s}` and to :math:`0.001~\text{m}^2/\text{s}`. Use a copy of the case for each value.

   - Before you run: calculate :math:`Re` and :math:`t_\text{visc}` for each value. Does the simulated duration of 0.5 s still make sense? At the higher viscosity, reduce :code:`deltaT` and save more frequently if the transient is too short to resolve. For :math:`\nu = 0.001~\text{m}^2/\text{s}`, increase :code:`endTime` in :code:`system/controlDict` (e.g. to 2 s).

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

   - Distinguish reaching steady flow in time from obtaining a result that changes little when the mesh is refined.

5. **A deep cavity: how many vortices?**

   - Make a copy of the base case and make the cavity four times as deep as it is wide (0.1 m wide, 0.4 m deep). In :code:`system/blockMeshDict`, change the y coordinates of the upper vertices. Which ones are they? Also change the number of cells in the y direction, so that the cells stay square.
   - Set :code:`endTime 1;` in :code:`system/controlDict`.
   - **Predict first:** before you run the case, sketch the flow you expect on a piece of paper. How many vortices will there be, and in which direction do they rotate?
   - Run the case and look at it in ParaView:

     - Color it by the velocity magnitude. What do you see?
     - Now switch to a logarithmic color scale: in the *Color Map Editor*, tick *Use Log Scale When Mapping Data To Colors* (click on the gear icon if you do not see this option), and set the range to 1e-7 to 1 with *Rescale to Custom Range*. What do you see now?
     - Make the flow paths visible with *Surface LIC*, a display mode that draws fine streaks along the flow direction:

       1. Apply a *Calculator* with the result name :code:`direction` and the expression :code:`U/max(mag(U), 1e-12)`. Away from zero velocity, this gives vectors of unit length, making weak vortices easier to see. The small lower bound avoids division by zero at stationary walls and stagnation points.
       2. Apply a second *Calculator* to the first one, with the result name :code:`log_speed` and the expression :code:`log10(max(mag(U), 1e-12))`. The lower bound avoids taking the logarithm of zero; the numerical values of :code:`U` are in m/s.
       3. In the toolbar, change the display mode from *Surface* to *Surface LIC*. In the *Properties* panel, set *Vectors* to :code:`direction`.
       4. Color by :code:`log_speed` and set the range to -6 to 0 with *Rescale to Custom Range*. The color now shows the speed: -1 means 0.1 m/s, -6 means 0.000001 m/s.

   - Questions to explore:

     - How many vortices are there? Was your sketch right?
     - The streaks show the local flow orientation, but not the direction of motion. In steady flow, streamlines also trace fluid paths. Think: the lid moves to the right, so in which direction does the top vortex rotate? Where two vortices touch, the fluid of both moves in the same direction, like two gears. So in which direction do the middle and the bottom vortex rotate?
     - Check your answer: use *Plot Over Line* along the vertical center line, from (0.05 0 0.005) to (0.05 0.4 0.005), and look at the horizontal velocity :code:`U_X` in the *SpreadSheet View*. Positive values mean flow to the right, negative values flow to the left. How much weaker is each vortex than the one above it?
     - What does this mean for the fluid at the bottom of a deep, narrow cavity, e.g. a deep fissure below a flowing current?
     - Optional: try a cavity that is two or three times as deep as it is wide. How does the number of vortices depend on the depth?

.. admonition:: Things to keep in mind

        - You may have to adjust the time step size to maintain a stable Courant number. 
        - When changing time step size, make sure to adjust the output interval accordingly.
        - Make a copy of the case for each variation, and show the results next to each other in ParaView (see "Comparing runs side by side" above).
        - Use the *Plot Over Line* filter in ParaView to extract velocity profiles and compare different cases.
        - Use the *Stream Tracer* filter in ParaView to visualize streamlines and vortex structure.




Expected Observations
---------------------

- **Baseline**: expect the center velocity to settle well before the characteristic diffusion time :math:`t_\text{visc} \sim 1` s. Test the approximate reference value of 0.07 s with your probe and 1 % criterion. The original output interval of 0.1 s misses much of this transient. A settled center velocity alone does not prove that the entire flow is steady.
- **High viscosity**: the transient is shorter and the Reynolds number is smaller. All three viscosity cases here are expected to remain laminar.
- **Low viscosity**: the transient is generally longer and inertial effects are stronger. Here :math:`Re = 100`; turbulence is not expected. Investigate growing oscillations for numerical causes before interpreting them as physical instabilities.
- **High lid velocity**: higher velocities and Reynolds number. The flow pattern may change, but the viscous diffusion time :math:`L^2/\nu` is unchanged; there is no general rule that the transient must become shorter or longer.
- **Deep cavity** (four times as deep as wide): three vortices on top of each other, with centers at about 0.38, 0.24 and 0.09 m height, i.e. about 1.4 cavity widths apart. Neighboring vortices rotate in opposite directions, like gears. The velocity decreases by a factor of several hundred between successive vortices: the largest horizontal velocities on the center line are about 0.19, 5e-4 and 1e-6 m/s. With a linear color scale, only the top vortex is visible; with *Surface LIC* on the direction field, all three show up as closed loops. Seen from the front, the top vortex rotates clockwise, the middle one counterclockwise, and the bottom one clockwise again. The fluid at the bottom is practically at rest. These are approximate reference values for comparison with your results; the weakest vortex is especially sensitive to mesh resolution and run duration. The alternating, rapidly weakening eddies are related to the viscous eddy sequences described by `Moffatt (1964) <https://www.damtp.cam.ac.uk/user/hkm2/PDFs/Moffatt_1964_JFM_18_1.pdf>`_, including the case of flow between parallel walls.
- **Fine mesh**: can capture finer details (e.g., secondary vortices) and reduce spatial discretization errors, at greater computational cost. Compare profiles to assess whether the solution is sufficiently resolved.

.. admonition:: Why is the center steady after about 0.07 s and not after 1 s?

    The estimate :math:`t_\text{visc} \sim L^2/\nu` gives a characteristic time scale, not the exact time at which a chosen velocity stops changing. A rough picture, not an exact calculation, helps to understand why the measured time is so much shorter:

    1. The motion of the lid does not have to diffuse across the whole cavity. The vortex is set up once the motion has reached about the center of the cavity, which is only half the distance, :math:`L/2`. Because the distance enters squared, half the distance needs only a quarter of the time.
    2. In the time :math:`t`, diffusion reaches about :math:`2\sqrt{\nu t}`, not :math:`\sqrt{\nu t}`. You know this from the cooling of the oceanic lithosphere, whose thickness grows as about :math:`2.3\sqrt{\kappa t}`.

    Together, the motion reaches the center when :math:`2\sqrt{\nu t} \approx L/2`, i.e. after

    .. math::

        t \approx \frac{L^2}{16\,\nu} = \frac{(0.1~\text{m})^2}{16 \cdot 0.01~\text{m}^2/\text{s}} \approx 0.06~\text{s},

    which is close to what we measure at the center. Diffusion has no sharp front, so the factor 1/16 is only approximate, and other parts of the cavity may need longer.

    At sufficiently low Reynolds number, increasing :math:`\nu` shortens the transient approximately in proportion to :math:`1/\nu`. Changing viscosity also changes :math:`Re`, so this proportionality need not be exact across all our runs. Use the scaling to predict trends, then measure the transient and check that your time step and output interval resolve it.
