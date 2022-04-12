thermo.correct();
rho=thermo.rho();
mu=thermo.mu();
Cp=thermo.Cp();
alphaP=thermo.alphaP();
betaT=thermo.betaT();
enthalpy=thermo.enthalpy();
// update heat flux
heatFlux = -kr*fvc::grad(T);

// pseudo gas generation
{
    forAll(alphaP,i)
    {
        if(T[i]>(200+273))
        {
            alphaP[i] = alphaP[i]*1.1;
        }       
    }
}

// update cp or rock according to temperature field
#include "updateCp_rock.H"
#include "update_permeability.H"

// Info<<cp_rock<<endl;