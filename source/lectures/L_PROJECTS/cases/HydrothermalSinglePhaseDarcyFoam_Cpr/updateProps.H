thermo.correct();
rho=thermo.rho();
mu=thermo.mu();
Cp=thermo.Cp();
alphaP=thermo.alphaP();
betaT=thermo.betaT();
enthalpy=thermo.enthalpy();
// update heat flux
heatFlux = -kr*fvc::grad(T);


// // 1. update internal value
// {
//     double h;//h: enthalpy
//     double T0,p0;
//     forAll(rho,i)
//     {
//         T0=T[i];p0=p[i];
//         if(p0>MAXP_FREESTEAM || p0<MINP_FREESTEAM)
//         {
//             Info<<"Error in freesteam: valid pressure range in freesteam is [1bar, 1000bar]. But pressure at the "<<i<<"th element is "<<p[i]<<endl;
//             exit(0);
//         }
//         if(T0>MAXT_FREESTEAM)
//         {
//             // Info<<"Error in freesteam: valid temperature range in freesteam is [5C, 800C]. But temperature at the "<<i<<"th element is "<<T[i]<<endl;
//             // Info<<"Debug: set low T to MAXT_FREESTEAM"<<endl;
//             T0=MAXT_FREESTEAM;
//             // exit(0);
//         }
//         if(T0<MINT_FREESTEAM)
//         {
//             // Info<<"Error in freesteam: valid temperature range in freesteam is [5C, 800C]. But temperature at the "<<i<<"th element is "<<T[i]<<endl;
//             // Info<<"Debug: set low T to MINT_FREESTEAM"<<endl;
//             T0=MINT_FREESTEAM;
//             // exit(0);
//         }
//         // Info<<"pressure: "<<p[i]<<"   Temperature: "<<T[i]<<endl;
//         SteamState S = freesteam_set_pT(p0, T0);
//         // rho[i]=freesteam_rho(S);
//         // mu[i]=freesteam_mu(S);
//         // cp_f[i]=freesteam_cp(S);
//         calculateProperties_h(S,h,rho[i],mu[i],betaT[i],alphaP[i],Cp[i]);
//         // calculateDlnRho_DlnT(p0,T0,dLnRho_dLnT[i]);
//     }
//     // 2. update boundary value
//     forAll(rho.boundaryField(), patchI)
//     {
//         forAll(rho.boundaryField()[patchI], faceI)
//         {
//             p0=p.boundaryField()[patchI][faceI];
//             T0=T.boundaryField()[patchI][faceI];
//             if(p0>MAXP_FREESTEAM || p0<MINP_FREESTEAM)
//             {
//                 Info<<"Error in freesteam: valid pressure range in freesteam is [1bar, 1000bar]. But pressure at the "<<faceI<<"th face of the "<<patchI<<"th patch is "<<p0<<endl;
//                 exit(0);
//             }
//             if(T0>MAXT_FREESTEAM)
//             {
//                 Info<<"Error in freesteam: valid temperature range in freesteam is [5C, 800C]. But temperature at the "<<faceI<<"th face of the "<<patchI<<"th patch is "<<T0<<endl;
//                 Info<<"Debug: set low T to MAXT_FREESTEAM"<<endl;
//                 T0=MAXT_FREESTEAM;
//                 // exit(0);
//             }
//             if(T0<MINT_FREESTEAM)
//             {
//                 Info<<"Error in freesteam: valid temperature range in freesteam is [5C, 800C]. But temperature at the "<<faceI<<"th face of the "<<patchI<<"th patch is "<<T0<<endl;
//                 Info<<"Debug: set low T to MINT_FREESTEAM"<<endl;
//                 T0=MINT_FREESTEAM;
//                 // exit(0);
//             }
//             // Info<<"boundary pressure: "<<p0<<"   Temperature: "<<T0<<endl;
            
//             SteamState S = freesteam_set_pT(p0, T0);
//             // rho.boundaryFieldRef()[patchI][faceI]=freesteam_rho(S);
//             // mu.boundaryFieldRef()[patchI][faceI]=freesteam_mu(S);
//             // cp_f.boundaryFieldRef()[patchI][faceI]=freesteam_cp(S);
//             calculateProperties_h(S,h,rho.boundaryFieldRef()[patchI][faceI],
//                                     mu.boundaryFieldRef()[patchI][faceI],
//                                     betaT.boundaryFieldRef()[patchI][faceI],
//                                     alphaP.boundaryFieldRef()[patchI][faceI],
//                                     Cp.boundaryFieldRef()[patchI][faceI]);
//             // calculateDlnRho_DlnT(p0,T0,dLnRho_dLnT.boundaryFieldRef()[patchI][faceI]);
//         }
//     }

    
// }
// // thermo.correct();   //update thermo
// // // Info<<"test thermo: "<<thermo.rho()<<endl;
// // // Info<<"write cp_f: "<<rho<<endl;
// // rho=thermo.rho();
// // cp_f=thermo.Cp();
// // mu=thermo.mu()*0.7;
// // beta_f=beta_f*300;
// // mu=mu*0.9;
