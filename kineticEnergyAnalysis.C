/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "kineticEnergyAnalysis.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kineticEnergyAnalysis, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
volScalarField Foam::kineticEnergyAnalysis::getTemporalKE()
{
    volScalarField temporalKE
    (
        IOobject
        (
            "temporalKE",
            U_.instance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (U_ & (U_ - U_.oldTime())/runTime_.deltaT())
    );
    temporalKE.write(runTime_.outputTime());


    scalar sumTemporalKE = gSum(temporalKE.primitiveFieldRef() * mesh_.V());
    Info << "sum of temporalDerivative KE: " << sumTemporalKE << endl;
    return temporalKE;
}

volScalarField Foam::kineticEnergyAnalysis::getConvectionKE()
{
    const volVectorField& C_ = mesh_.lookupObject<volVectorField>(convName_);


    volScalarField convKE
    (
        IOobject
        (
            "convKE",
            U_.instance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (U_ & C_)
    );
    convKE.write(runTime_.outputTime());


    scalar sumConvKE = gSum(convKE.primitiveFieldRef() * mesh_.V());
    Info << "sum of convection KE: " << sumConvKE << endl;
    return convKE;
}

volScalarField Foam::kineticEnergyAnalysis::getPGradKE()
{
    // Approach 3: (use least square pressure gradient)
    // volVectorField pGradLS_
    //               ("pGradLS", mesh_.lookupObject<volVectorField>("pGradLS"));

    volScalarField pGradKE
    (
        IOobject
        (
            "pGradKE",
            U_.instance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (U_ & fvc::grad(p_)) // Approach 1
        // (U_ & fvc::reconstruct(fvc::snGrad(p_)* mesh_.magSf())) // Approach 2
        // (U_ & pGradLS_) // Approach 3
    );
    pGradKE.write(runTime_.outputTime());

    scalar sumPGradKE = gSum(pGradKE.primitiveFieldRef() * mesh_.V());
    Info << "sum of pGrad KE: " << sumPGradKE << endl;
    return pGradKE;
}


volScalarField Foam::kineticEnergyAnalysis::getDissipationKE()
{
    // Calculate dissipation
    tmp<volScalarField> nuEff_
    (
        mesh_.lookupObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        ).nuEff()
    );

    volScalarField dissipKE
    (
        IOobject
        (
            "dissipKE",
            U_.instance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_ & (  fvc::laplacian(nuEff_.ref(), U_)
              + fvc::div(nuEff_.ref()*dev(T(fvc::grad(U_))))
             )
    );
    dissipKE.write(runTime_.outputTime());
    nuEff_.clear();

    scalar sumDissipKE = gSum(dissipKE.primitiveFieldRef() * mesh_.V());
    Info << "sum of dissipation KE: " << sumDissipKE << endl;
    return dissipKE;
}


void Foam::kineticEnergyAnalysis::analyzeKEBalance()
{
    // volScalarField temporalKE
    // (
    //     IOobject
    //     (
    //         "temporalKE",
    //         U_.instance(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     (U_ & (U_ - U_.oldTime())/runTime_.deltaT())
    // );
    // temporalKE.write(runTime_.outputTime());


    // scalar sumTemporalKE = gSum(temporalKE.primitiveFieldRef() * mesh_.V());
    // Info << "sum of temporalDerivative KE: " << sumTemporalKE << endl;



    // volScalarField convKE
    // (
    //     IOobject
    //     (
    //         "convKE",
    //         U_.instance(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     (U_ & C_)
    // );
    // convKE.write(runTime_.outputTime());


    // scalar sumConvKE = gSum(convKE.primitiveFieldRef() * mesh_.V());
    // Info << "sum of convection KE: " << sumConvKE << endl;


    // volScalarField p_ ("p", mesh_.lookupObject<volScalarField>("p") );

    // // Approach 3: (use least square pressure gradient)
    // // volVectorField pGradLS_("pGradLS", mesh_.lookupObject<volVectorField>("pGradLS"));

    // volScalarField pGradKE
    // (
    //     IOobject
    //     (
    //         "pGradKE",
    //         U_.instance(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     (U_ & fvc::grad(p_)) // Approach 1
    //     // (U_ & fvc::reconstruct(fvc::snGrad(p_)* mesh_.magSf())) // Approach 2
    //     // (U_ & pGradLS_) // Approach 3
    // );
    // pGradKE.write(runTime_.outputTime());

    // scalar sumPGradKE = gSum(pGradKE.primitiveFieldRef() * mesh_.V());
    // Info << "sum of pGrad KE: " << sumPGradKE << endl;




    // volScalarField dissipKE
    // (
    //     IOobject
    //     (
    //         "dissipKE",
    //         U_.instance(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     U_ & (  fvc::laplacian(nuEff_.ref(), U_)
    //           + fvc::div(nuEff_.ref()*dev(T(fvc::grad(U_))))
    //          )
    // );
    // dissipKE.write(runTime_.outputTime());
    // nuEff_.clear();

    // scalar sumDissipKE = gSum(dissipKE.primitiveFieldRef() * mesh_.V());
    // Info << "sum of dissipation KE: " << sumDissipKE << endl;



    volScalarField KEbalance
    (
        IOobject
        (
            "KEbalance",
            U_.instance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0)
    );

    KEbalance.primitiveFieldRef() =
        (  getTemporalKE()
         + getConvectionKE()
         + getPGradKE()
         - getDissipationKE()
         );
    KEbalance.write(runTime_.outputTime());

    scalar sumKEbalance = gSum(KEbalance.primitiveFieldRef() * mesh_.V());
    Info << "sum of KE balance: " << sumKEbalance << endl;

    return;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticEnergyAnalysis::kineticEnergyAnalysis
(
    const volVectorField& U,
    const volScalarField& p,
    const word convName
)
:
    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    // Set the pointer to the velocity, flux and pressure-correction field
    U_(U),
    p_(p),
    convName_(convName)
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticEnergyAnalysis::~kineticEnergyAnalysis()
{
    Info << "kineticEnergyAnalysis Destructor" << endl;
}


// ************************************************************************* //
