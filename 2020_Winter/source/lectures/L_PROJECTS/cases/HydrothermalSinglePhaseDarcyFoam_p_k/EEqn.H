/**
 * @file EEqn.H
 * @brief Solve energy conservation equation.
 * 
 * - **Energy conservation equation in terms of temperature**
 * \f[ \left( {\phi {\rho _f}{c_{pf}} + \left( {1 - \phi } \right){\rho _r}{c_{pr}}} \right)\frac{{\partial T}}{{\partial t}} = \nabla  \cdot \left( {{\lambda _r}\nabla T} \right) - {\rho _f}{c_{pf}} \vec{v} \cdot \nabla T + \frac{{{\mu _f}}}{k}{\vec{v}^2} - {\left( {\frac{{\partial \;ln\rho }}{{\partial \;lnT}}} \right)_p}\frac{{Dp}}{{Dt}} \f]
 * where \f$ \left( {\frac{{\partial \;ln\rho }}{{\partial \;lnT}}} \right)_p = -T\alpha_p \f$, 
 * \f$ \alpha_p \f$ is the fluid thermal expansivity.
 * 
 * - **Energy conservation equation in terms of enthalpy**
 * 
 * @version 1.0
 * @author Zhikui Guo (zguo@geomar.de)
 * @date 2019-10-14
 * 
 * @copyright Copyright (c) 2019
 * 
 */

if (thermo.he().name()=="T")
{
    fvScalarMatrix TEqn
    (
        (porosity*rho*Cp+(1.0-porosity)*rho_rock*cp_rock)*fvm::ddt(T)
        +fvm::div(phi*fvc::interpolate(Cp),T)
        == 
        fvm::laplacian(kr,T)
        +
        fvm::Sp(fvc::div(phi*fvc::interpolate(Cp)),T)
        +
        mu/permeability*magSqr(U)
        +
        fvm::Sp(alphaP*(porosity*fvc::ddt(p)+(U & fvc::grad(p))),T)
    );

    TEqn.solve();
    // T.correctBoundaryConditions();
}else
{
    FatalErrorInFunction
            << "The energy equation in terms of enthalpy is still under developing\n" 
            << "Please use keyword of temperture for energy entry in thermophysicalProperties."
            << abort(FatalError);
}