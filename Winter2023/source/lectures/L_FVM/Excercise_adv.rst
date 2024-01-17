.. include:: /include.rst_

.. _Excercise_adv:

==============================
Excercise on advection schemes
==============================


In this excercise we will compare the performance of the different advection schemes. We will use the general advection-diffusion from previsou sections:

.. math::
    \frac{\partial T}{\partial t} = \nabla \cdot \left( D \nabla T \right) - \nabla \cdot \left( \vec{U} T \right)


Solver
-------
OpenFoam has a simplified solver for the advection-diffusion equation called :code:`scalarTransportFoam`, which uses a prescribed velocity field to move the scalar field around.

As always, a good starting point is find a suitable tutorial, copy it over to your own case and start modifying it. In this case, we will use the tutorial :code:`basic/scalarTransportFoam/pitzDaily` as a starting point. Copy the case to your own case directory and give it a reasonable name.

.. code-block:: bash

    cp -r $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily $HOME/HydrothermalFoam_runs


The solver is set up to solve the advection-diffusion equation for a scalar field called :code:`T`. The velocity field is prescribed in the file :code:`0/U` and the diffusion coefficient is set in :code:`constant/transportProperties`. The initial condition is set in :code:`0/T`.

Mesh and initial conditions
------------------

Mesh
^^^^
We will use a 2D box with dimensions :math:`-0.5 < x < 0.5` and :math:`-0.5 < y < 0.5`. You will need to modify the blockMeshDict file to set up the mesh. The mesh should be uniform with 100 cells in each direction. The mesh should be centered around the origin.

Temperature
^^^^^^^^^^^^

The scalar field initial conditions are set in the file :code:`0/T`. We will use a Gaussian distribution with a standard deviation of :math:`\sigma = 0.1` and a maximum value of :math:`T_0 = 2`. It should be centered at x. Let's put the intitial gaussian at :math:`x_0 = 0` and :math:`y_0 = 0.25`.

.. math::
    T(x,y) = T_0 \exp \left( - \frac{(x-x_0)^2 + (y-y_0)^2}{\sigma^2} \right)    


Velocity
^^^^^^^^^

We will fist look into solid body rotation, in which the initial gaussian will be rotated in clockwise direction without any shearing. The respective velocity field is given by:

.. math::

    Vx(x,y)=y 
    
.. math::
    Vy(x,y)=-x 


Implementation
^^^^^^^^^^^^^^
To implement the initial conditions and prescribed velocity field, you will need to modify the files :code:`0/T` and :code:`0/U`. We will again use codestream statements for this.


.. code-block:: foam
    :linenos:
    :emphasize-lines: 36-61

    /*--------------------------------*- C++ -*----------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     | Website:  https://openfoam.org
      \\  /    A nd           | Version:  7
       \\/     M anipulation  |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      T;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 0 0 1 0 0 0];

    //internalField   uniform 278.15;
    internalField #codeStream
    {
        codeInclude
        #{
            #include "fvCFD.H"
        #};
        codeOptions
        #{
            -I$(LIB_SRC)/finiteVolume/lnInclude \
            -I$(LIB_SRC)/meshTools/lnInclude
        #};
        codeLibs 
        #{
            -lmeshTools \
            -lfiniteVolume 
        #};
        localCode
        #{
            static double calTemperature(const scalar sigma, const scalar maxT, scalar x, scalar y)
            {
                return ???;

            }
        #};
        code 
        #{
            const IOdictionary& d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());
            scalarField T(mesh.nCells(), 0);
            scalar sigma = .1, maxT = 2;
    
            forAll(T, i)
            {
                const scalar x = mesh.C()[i][0];
                const scalar y = mesh.C()[i][1];
                T[i]=calTemperature(sigma, maxT,x, y);            
            }
        

            writeEntry(os, "", T); //
        #};
    };



    boundaryField
    {
        left
        {
            type            zeroGradient;
        }

        right
        {
            type            zeroGradient;
        }

        top
        {
            type            zeroGradient;
        }

        bottom
        {
            type            zeroGradient;
        }

        frontAndBack
        {
            type            empty;
        }
    }

    // ************************************************************************* //


We can do something similar for the velocity field. The velocity field is set in the file :code:`0/U`. Again, we will use codestream statements to set the velocity field.

.. code-block:: foam
    :linenos:
    :emphasize-lines: 36-55

    /*--------------------------------*- C++ -*----------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     | Website:  https://openfoam.org
      \\  /    A nd           | Version:  7
        \\/     M anipulation  |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volVectorField;
        location    "0";
        object      U;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 1 -1 0 0 0 0];

    internalField #codeStream
    {
        codeInclude
        #{
            #include "fvCFD.H"
        #};
        codeOptions
        #{
            -I$(LIB_SRC)/finiteVolume/lnInclude \
            -I$(LIB_SRC)/meshTools/lnInclude
        #};
        codeLibs 
        #{
            -lmeshTools \
            -lfiniteVolume 
        #};
        code 
        #{
            const double pi = 3.141592653589793;
            const IOdictionary& d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());
            volVectorField& U = const_cast<volVectorField&>(
                mesh.lookupObject<volVectorField>("U")
            );

            forAll(U, celli)
            {
                const point& coord = mesh.C()[celli];
                U[celli].x() = ???;
                U[celli].y() = ???;
                U[celli].z() = 0;
            }
            
            writeEntry(os, "", U); //
        #};
    };


    boundaryField
    {
        left
        {
            type            zeroGradient;
        }
        right
        {
            type            zeroGradient;
        }
        top
        {
            type            zeroGradient;
        }
        bottom
        {
            type            zeroGradient;
        }
        frontAndBack
        {
            type            empty;
        }
    }


    // ************************************************************************* //


Setting up the case
^^^^^^^^^^^^^^^^^^^

Since we want to look at advection without any physical diffusion, we also need to set the diffusion coefficient to something small. This happens in :code:`constant/transportProperties`:


.. code-block:: foam
    :linenos:
    :emphasize-lines: 36-61

    /*--------------------------------*- C++ -*----------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     | Website:  https://openfoam.org
      \\  /    A nd           | Version:  7
        \\/     M anipulation  |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      transportProperties;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    DT              DT [0 2 -1 0 0 0 0] 1e-9;


    // ************************************************************************* //

Finally, we need to set the time step and simulation time. This happens in :code:`system/controlDict`.

.. code-block:: foam
    :linenos:
    :emphasize-lines: 26-32


    /*--------------------------------*- C++ -*----------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     | Website:  https://openfoam.org
      \\  /    A nd           | Version:  7
        \\/     M anipulation  |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      controlDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    application     scalarTransportFoam;

    startFrom       startTime;

    startTime       0;

    stopAt          endTime;

    endTime         20;

    deltaT          0.0025;

    writeControl    timeStep;

    writeInterval   20;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;


    // ************************************************************************* //


We here use a constant rund time of 20 and a time step length of :math:`\Delta t = 0.0025`; we save the output every 20 steps. To run the case, call :code:`blockMesh` and then the solver :code:`scalarTransportFoam`. 

Excercise
----------

Finally, change the advection scheme (e.g. vanLeer, upwind, MUSCL) and compare the results in paraview! A good way is to make multiple directions for the individual cases.

Now, change the velocity field to a shear shell and repeat the exercise. 

.. math::
    
    Vx(x,y) = -sin(\pi*(X+0.5))*cos(\pi*(Y+0.5))
    
.. math:: 
    
    Vy(x,y) = cos(\pi*(X+0.5))*sin(\pi*(Y+0.5))


