/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "functionObjectTemplate.H"
#include "fvCFD.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(plumeTemperatureFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    plumeTemperatureFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 458a6d83b7c6a4a9daeed2389c4cb96ebf80410e
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void plumeTemperature_458a6d83b7c6a4a9daeed2389c4cb96ebf80410e(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const fvMesh& plumeTemperatureFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plumeTemperatureFunctionObject::plumeTemperatureFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

plumeTemperatureFunctionObject::~plumeTemperatureFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool plumeTemperatureFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read plumeTemperature sha1: 458a6d83b7c6a4a9daeed2389c4cb96ebf80410e\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool plumeTemperatureFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute plumeTemperature sha1: 458a6d83b7c6a4a9daeed2389c4cb96ebf80410e\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool plumeTemperatureFunctionObject::write()
{
    if (false)
    {
        Info<<"write plumeTemperature sha1: 458a6d83b7c6a4a9daeed2389c4cb96ebf80410e\n";
    }

//{{{ begin code
    #line 54 "/Users/lruepke/Research/Lectures/HTF_lecture/Summer2022/source/lectures/L05/cases/Driesner2010/HS86kw/5E-14/system/controlDict.functions.calPlumeT"
//get maximum tempeature on the top boundary
            label patchID = mesh().boundaryMesh().findPatchID("top"); 
            const volScalarField& T = mesh().lookupObject<volScalarField>("T");
            // write vent temperature
            std::ofstream fout("ventT.txt",std::ofstream::app);
            // Info<<"Plume Temperature: "<<mesh().time().value()/31536000<<"\t"<<Foam::gMax(T.boundaryField()[patchID])-273.15<<" C"<<endl;
            fout<<mesh().time().value()/31536000<<"\t"<<Foam::gMax(T.boundaryField()[patchID])-273.15<<std::endl;
            fout.close();
//}}} end code

    return true;
}


bool plumeTemperatureFunctionObject::end()
{
    if (false)
    {
        Info<<"end plumeTemperature sha1: 458a6d83b7c6a4a9daeed2389c4cb96ebf80410e\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

