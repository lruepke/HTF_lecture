.. include:: /include.rst_

.. _Advection1:

====================
The advection term
====================

In the previous sessions, we have explored in details how openFoam solves the diffusion equation. In this session, we will explore how openFoam solves the advection equation.

Do you remember the governing equations of hydrothermal convection?

.. math::
    :label: eq:pressure_adv_fvm

    \varepsilon \rho_f \left( \beta_f \frac{\partial p}{\partial t} - \alpha_f \frac{\partial T}{\partial t} \right) = \nabla \cdot \left( \rho_f \frac{k}{\mu_f} (\nabla p - \rho_f \vec{g}) \right)


.. math::
    :label: eq:temperature_adv_fvm

    (\varepsilon \rho_f C_{pf} + (1-\varepsilon)\rho_r C_{pr})\textcolor{green}{\frac{\partial T}{\partial t}} = \textcolor{orange}{\nabla \cdot (\lambda_r \nabla T)} \textcolor{purple}{- \rho_f C_{pf} \vec{U}\cdot \nabla T} + \frac{\mu_f}{k} \parallel \vec{U} \parallel ^2 - \left( \frac{\partial ln \rho_f}{\partial ln T} \right)_p \frac{Dp}{Dt}


We have already looked into the numerical FVM implementation of the time derivative (green term)and the diffusion term (orange term). The purple term is the transport, so-called advection term. It handles the transport of variable by the flowing fluid. In the hydrothermal convection case, we have transported the temperature, :math:`T`, with the fluid.

In a more general form, a transient advection diffusion equation can be written as:

.. math::
    \frac{\partial T}{\partial t} = \nabla \cdot \left( D \nabla T \right) - \nabla \cdot \left( \vec{U} T \right)

where :math:`T` is the transported scalar, :math:`D` is the diffusion coefficient and :math:`\vec{U}` is the velocity field. 

.. admonition:: Forms of the advection term
    
    Did you notice the difference between the temperature equation for hydrothermal convection and the general advection equation? In the form given for hydrothermal convection, the specific heat :math:`c_p`, the density :math:`\rho`, and the velocity :math:`\vec{u}` are outside the divergence, while in the general form they are inside (due to the conservative form of the equation). The advection term is also called the convective term.

    The reason is not that these variables are constant but lies in the way the energy conservation equation is derived. If it is written in terms of energy (e.g. enthalpy or internal energy), we would have the respective terms included in the divergence. If we use the temperature form, using the specific heat, those terms end up outside the divergence. You can check :cite:`Guo2020` for details.


Resolving advection
-------------------

The core of the advection or transport term solving process involves the discretization of the convective fluxes across the control volume faces. This discretization is crucial as it directly influences the stability and accuracy of the simulation. OpenFOAM provides a range of discretization schemes, such as the upwind, linear, and higher-order schemes like QUICK and cubic, which are discussed and presented in :cite:`moukalled2016finite`. The choice of scheme can significantly affect the numerical diffusion and the overall solution's fidelity.

For instance, the upwind scheme, while robust and simple, can introduce considerable numerical diffusion, leading to a less accurate representation of sharp gradients. On the other hand, higher-order schemes offer improved accuracy but can lead to numerical instabilities if not handled properly. You can find an introduction to different advection scheme in the companian lecture on finite elements and general numerical methods (`Link <https://lruepke.github.io/HTF_lecture/winter2022/lecture1/fdm_advection.html>`_)

Moreover, OpenFOAM incorporates the concept of boundedness through the use of Total Variation Diminishing (TVD) schemes to ensure that the physical and numerical bounds are respected, especially in cases involving steep gradients or when dealing with scalar quantities like concentration or temperature.

In essence, the numerics of solving the advection/transport term in OpenFOAM involve a careful balance between accuracy, stability, and computational efficiency, with the choice of discretization schemes and iterative solvers playing pivotal roles in the outcome of the CFD simulations.

Implementation in OpenFoam
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The way the advection term is solved is again handled in the file :code:`system/fvSchemes`.

.. code-block::
      :linenos:
      :emphasize-lines: 28-32,34-38
      
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
        object      fvSchemes;
      }
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

      ddtSchemes
      {
        default         Euler;
      }

      gradSchemes
      {
        default         Gauss linear;
      }

      divSchemes
      {
         default         none;
         div(phi,T)      Gauss MUSCL grad(T);
      }

      laplacianSchemes
      {
        default         none;
        laplacian(DT,T) Gauss linear corrected;
      }

       interpolationSchemes
      {
        default         linear;
      }

      snGradSchemes
      {
        default         corrected;
      }


      // ************************************************************************* //


.. tip::
    There are many different advection schemes implemented in OpenFoam. You can get a list of all available schemes by changing the advection scheme to something invalid like :code:`abc`. Openfoam will throw an error and list all available schemes.

