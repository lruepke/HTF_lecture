/**
 * @file initpEqn.H
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Calculate initial hydrostatic pressure.
 * 
 * \f[ \nabla  \cdot \left( {{\rho _f}\frac{k}{{{\mu _f}}}\left( {\nabla p - {\rho _f}\vec{g}} \right)} \right)=0 \f]
 * @version 1.0
 * @date 2019-10-14
 * 
 * @copyright Copyright (c) 2019
 * 
 */
{
    volScalarField rAU(permeability/mu); //K/mu
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));//rho/A on surface
    volVectorField HbyA(U*0);

    surfaceScalarField phig("phig",(fvc::interpolate(rho)*rhorAUf * g) & mesh.Sf());
    
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        phig
    );
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf);

    // fvScalarMatrix p_rghDDtEqn
    // (
    //     porosity*rho*beta_f*fvm::ddt(p)
    //     -porosity*rho*alpha_f*fvc::ddt(T)
    //     +
    //     fvc::div(phiHbyA)
    // );
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA) - fvm::laplacian(rhorAUf, p)
        );
        pEqn.solve();
        if (pimple.finalNonOrthogonalIter())
        {
            // option 1: calculate velocity using Darcy's law directly
            // U=-permeability/mu*(fvc::grad(p)-rho*g);

            // option 2: using flux reconstruct velocity
            phi = phiHbyA + pEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + pEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
        }
    }
}
