/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------*/
 
/**
 * @file HydrothermalSinglePhaseDarcyFoam.C
 * @brief The main program of HydrothermalSinglePhaseDarcyFoam solver.
 * 
 * which is developed to modeling **porous flow** and **heat transfer** for solving 
 * hydrothermal circulation problem (e.g. [Hasenclever et al., 2014](https://doi.org/10.1038/nature13174)).
 * 
 * \dotfile HydrothermalSinglePhaseDarcyFoam.dot
 * @author Zhikui Guo (zguo@geomar.de)
 * @version 1.0 
 * @date 2019-10-14
 * @copyright Copyright (c) 2020 by [Zhikui Guo](https://www.modernfig.cn) and [Lars Rüpke](https://www.geomar.de/index.php?id=lruepke)
 * 
 */

#include "fvCFD.H"
#include "pimpleControl.H" 
#include "simpleMatrix.H"  //for investigating coefficients matrix
// new defined thermophysical model: hydroThermo, htHydroThermo, see libraries
#include "hydroThermo.H"
//---------------------------------------------
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    // initialize pressure: correct initial hydrostatic pressure is very important!!!
    if(runTime.timeName()=="0")
    {
        Info<<"Initialize hydrostatic pressure"<<endl;
        #include "initpEqn.H"
        #include "updateProps.H"
        #include "initpEqn.H"
        // #include "updateProps.H"
        Info<<"Initialize hydrostatic pressure end"<<endl;
    }
    // Get some constant values
    const scalar endTime( readScalar(runTime.controlDict().lookup("endTime")));
    const scalar sec_year = 86400*365;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        
        #include "EEqn.H"
        #include "pEqn.H"
        #include "updateProps.H"
        
        // ---------------------------------------

        runTime.write();
        
        #include "reportTotalExecutionTime.H"
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
