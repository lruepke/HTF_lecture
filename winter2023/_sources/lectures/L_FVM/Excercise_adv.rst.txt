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

Initial conditions
------------------

We will use a 2D bos with dimensions :math:`-0.5 < x < 0.5` and :math:`-0.5 < y < 0.5`. You will need to modify the blockMeshDict file to set up the mesh. The mesh should be uniform with 100 cells in each direction. The mesh should be centered around the origin.

The scalar field initial conditions are set in the file :code:`0/T`. We will use a Gaussian distribution with a standard deviation of :math:`\sigma = 0.1` and a maximum value of :math:`T_0 = 2`. It should be centered at x. 

.. math::
    T(x,y) = T_0 \exp \left( - \frac{(x-x_0)^2 + (y-y0)^2}{\sigma^2} \right)    



We will look into two different velocity fields:

.. math::

    Vx(x,y)=y 
    
.. math::
    Vy(x,y)=-x 

and

.. math::
    
    Vx(x,y) = -sin(\pi*(X+0.5))*cos(\pi*(Y+0.5))
    
.. math:: 
    
    Vy(x,y) = cos(\pi*(X+0.5))*sin(\pi*(Y+0.5))



To implmenent the intitial conditions and prescribed velocity field, you will need to modify the files :code:`0/T` and :code:`0/U`. We will again use codestream statements for this.


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


Try to implement somthing similar for the velocity and test different advection schemes. 

