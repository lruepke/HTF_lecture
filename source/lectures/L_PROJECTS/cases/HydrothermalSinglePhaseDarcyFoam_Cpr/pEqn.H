/**
 * @file pEqn.H
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief 
 * 
 * \f[ 
 * \phi {\rho _f}\left( {{\beta _f}\frac{{\partial p}}{{\partial t}} - {\alpha _f}\frac{{\partial T}}{{\partial t}}} \right) = \nabla  \cdot \left( {{\rho _f}\frac{k}{{{\mu _f}}}\left( {\nabla p - {\rho _f}\vec g} \right)} \right)
 * \f]
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

    fvScalarMatrix p_rghDDtEqn
    (
        porosity*rho*betaT*fvm::ddt(p)
        -porosity*rho*alphaP*fvc::ddt(T)
        +fvc::div(phiHbyA)
    );
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            p_rghDDtEqn - fvm::laplacian(rhorAUf, p)
        );
        pEqn.solve();
        if (pimple.finalNonOrthogonalIter())
        {
            // option 2: using flux reconstruct velocity, magic function of reconstruct
            phi = phiHbyA + pEqn.flux();
            U = HbyA + rAU*fvc::reconstruct((phig + pEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
        }
    }
}