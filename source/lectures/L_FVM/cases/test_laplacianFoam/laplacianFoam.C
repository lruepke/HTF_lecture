/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"

    // 1. read mesh
    #include "createMesh.H"
    #include "createFields.H"
    // access internal faces
    surfaceVectorField Cf = mesh.Cf();
    surfaceVectorField Sf = mesh.Sf();
    surfaceScalarField S = mesh.magSf();
    Info<<"====Cf start"<<endl;
    forAll(Cf, iface)
    {
        Info<<iface<<": face center "<<Cf[iface]<<endl;
        Info<<iface<<": face area vector "<<Sf[iface]<<endl;
        Info<<iface<<": face area "<<S[iface]<<endl;
        Info<<iface<<": face delta coeff "<<mesh.deltaCoeffs()[iface]<<endl;
        Info<<iface<<": coeff(D*magSf*deltacoeff) "<<mesh.deltaCoeffs()[iface]*DT*S[iface]<<endl;
    }
    // faceCentres();
    // Info<<"====Cf end"<<endl;

    // 2. access boundary mesh
    // const fvBoundaryMesh& boundaryMesh = mesh.boundary(); 
    // forAll(boundaryMesh, patchI)
    // {
    //     const fvPatch& patch = boundaryMesh[patchI];
    //     forAll(patch, faceI)
    //     {
    //         vector faceNormal = patch.Sf()[faceI]; 
    //         scalar faceArea = patch.magSf()[faceI]; 
    //         vector unitFaceNormal = patch.nf()()[faceI]; 
    //         vector faceCenter = patch.Cf()[faceI]; 
    //         label owner = patch.faceCells()[faceI];
    //     } 
    // }
    
    // 3. access boundary field, boundary field coefficients, 
    forAll(T.boundaryField(), patchI)
    {
        Info<<"Boundary patch: "<<mesh.boundary()[patchI].name()<<endl;
        Info<<"Is coupled ? "<<mesh.boundary()[patchI].coupled()<<endl;
        Info<<"gradientInternalCoeffs of field T "<<endl;
        Info<<T.boundaryField()[patchI].gradientInternalCoeffs()<<endl; //Diagonal coeff [A]
        Info<<"gradientBoundaryCoeffs of field T "<<endl;
        Info<<T.boundaryField()[patchI].gradientBoundaryCoeffs()<<"\n"<<endl; //source coeff, [B]
        forAll(T.boundaryField()[patchI], faceI)
        {
            // fvm.internalCoeffs()[patchi] = pGamma*pvf.gradientInternalCoeffs();
            // T.boundaryField()[patchI][faceI].gradientInternalCoeffs();
        }
    }

    simpleControl simple(mesh);

    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        fvScalarMatrix Laplacian(fvm::laplacian(DT, T));
        // Info<<"fvm::laplacian(DT, T): "<<"\n"
        //     <<"\tLower"<<Laplacian.lower()<<"\n"
        //     <<"\tDiagonal"<<Laplacian.diag()<<"\n"
        //     <<"\tUpper"<<Laplacian.upper()<<"\n"
        //     <<"\tinternalCoeffs"<<Laplacian.internalCoeffs()<<"\n"
        //     <<"\tboundaryCoeffs"<<Laplacian.boundaryCoeffs()<<"\n"
        //     <<"\tSource"<<Laplacian.source()<<"\n"
        //     <<endl;
        // Note that the boundary conditions related source is not assembled into the Laplacian after calling  fvm::laplacian(DT, T)
        // assembly of BCs related source will happed in the .solve() function, in the solve function, the addBoundarySource function will be called to do this.
        // so if you print the Laplacian.source(), you will find it is a zero array.
        // Laplacian.addBoundarySource(Laplacian.source(), false); // addBoundarySource is a protected function in fvMatrix Class
        Info<<"fvm::laplacian: "<<Laplacian<<endl;
        
        while (simple.correctNonOrthogonal())
        {
            // Info<<"ddt"<<endl;
            // Info<<fvm::ddt(T)<<endl;
            // Info<<"laplacian"<<endl;
            // Info<<fvm::laplacian(DT, T)<<endl;
            // Info<<"fvOptions(T)"<<endl;
            // Info<<fvOptions(T)<<endl;
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
            //  ==
            //     fvOptions(T)
            );

            // Info<<"TEqn source: \n"<<TEqn<<endl;
            // Info<<"TEqn"<<endl;
            // Info<<TEqn<<endl;
            // fvOptions.constrain(TEqn);
            TEqn.solve();
            // fvOptions.correct(T);
        }

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
