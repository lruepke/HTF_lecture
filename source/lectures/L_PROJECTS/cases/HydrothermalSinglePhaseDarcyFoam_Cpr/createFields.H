autoPtr<hydroThermo> pThermo(hydroThermo::New(mesh));
hydroThermo& thermo = pThermo();
// new defined thermophysical model: htHydroThermo, supporting temperature and enthalpy for energy equation
thermo.validate(args.executable(), "h","T"); 

// --------------non const access variables------------------------------------------
// pressure: only solved from pressure equation
volScalarField& p=thermo.p();
// Temperature field: solved from diffusion equation and updated in ThermalEquilibrium function
volScalarField& T=thermo.T();
// ----------------------------------------------------------------------------------

// Darcy Velocity
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U",dimensionSet(0,1,-1,0,0,0,0),vector::zero)
);

// Permeability
volScalarField permeability
(
    IOobject
    (
        "permeability",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
    
// change log: rock properties are moved from transportProperties to porousMedia sub-dict of thermophysicalProperties
dimensionedScalar rho_rock=thermo.rho_rock();
dimensionedScalar cp0_rock = thermo.cp_rock();
dimensionedScalar porosity0 = thermo.porosity();
dimensionedScalar kr = thermo.kr();

volVectorField heatFlux
(
    IOobject
    (
        "heatFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -kr*fvc::grad(T)
);

// Porosity, one can specify uniform porosity in thermophysicalProperties, and also specify it like permeability
volScalarField porosity
(
    IOobject
    (
        "porosity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    porosity0
);
// Cp_rock, one can specify uniform value in thermophysicalProperties, and also specify it like permeability in 0 folder
volScalarField cp_rock
(
    IOobject
    (
        "cp.rock",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    cp0_rock
);

// Fluid density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);
    
// fluid dynamic viscosity
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.mu()
);

//  specific heat of fluid
volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()
);
//  specific enthalpy
volScalarField enthalpy
(
    IOobject
    (
        "enthalpy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.enthalpy()
);

//  thermal expansivity
volScalarField alphaP
(
    IOobject
    (
        "alphaP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.alphaP()
);

//  compresibility
volScalarField betaT
(
    IOobject
    (
        "betaT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.betaT()
);

#include "readGravitationalAcceleration.H"

#include "compressibleCreatePhi.H"
