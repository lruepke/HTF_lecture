
{
    // volScalarField rAU(permeability/mu); //K/mu
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_rock));//rho on surface
    volVectorField HbyA(U*0);

    surfaceScalarField phig("phig",(fvc::interpolate(rho_rock)*rhorAUf * g) & mesh.Sf());
    
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        phig
    );
    // Update the pressure BCs to ensure flux consistency
    // constrainPressure(p, rho_rock, U, phiHbyA, rhorAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix plEqn
        (
            fvc::div(phiHbyA) - fvm::laplacian(rhorAUf, pl)
        );
        plEqn.solve();
    }
}
