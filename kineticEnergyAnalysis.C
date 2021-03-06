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

tmp<volScalarField> Foam::kineticEnergyAnalysis::getTemporalKE()
{
    tmp<volScalarField> tTemporalKE
    (
        new volScalarField
        (
            IOobject
            (
                "temporalKE",
                U_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            // (U_ & fvc::ddt(U_))
            // Inaccurate: U^(n+1) * (U^(n+1) - U^n)/dt
            // (U_ & (U_ - U_.oldTime()) /runTime_.deltaT() )

            // Accurate: ((U^(n+1))^2 - (U^n)^2)/dt
            0.5*((magSqr(U_) - magSqr(U_.oldTime())) /runTime_.deltaT())
        )
    );
    volScalarField& temporalKE = tTemporalKE.ref();
    temporalKE.write(runTime_.outputTime());

    // volume weighted average
    avgDdtKE_ = temporalKE.weightedAverage(mesh_.V());
    Info << "\tavg of temporalDerivative KE: " << avgDdtKE_.value()<< endl;

    return tTemporalKE;
}

tmp<volScalarField> Foam::kineticEnergyAnalysis::getConvectionKE()
{
    tmp<volScalarField> tConvKE
    (
        new volScalarField
        (
            IOobject
            (
                "convKE",
                U_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (U_ & fvc::div(phi_,U_))
        )
    );
    volScalarField& convKE = tConvKE.ref();
    convKE.write(runTime_.outputTime());

    // volume weighted average
    avgConvKE_ = convKE.weightedAverage(mesh_.V());
    Info << "\tavg of convection KE: " << avgConvKE_.value() << endl;

    return tConvKE;
}

tmp<volScalarField> Foam::kineticEnergyAnalysis::getPGradKE()
{
    // Approach 3: (use least square pressure gradient)
    // volVectorField pGradLS_
    //               ("pGradLS", mesh_.lookupObject<volVectorField>("pGradLS"));

    tmp<volScalarField> tPGradKE
    (
        new volScalarField
        (
            IOobject
            (
                "pGradKE",
                U_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            // mesh_,
            // dimensionedScalar("", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0)
            -(U_ & fvc::grad(p_)) // Approach 1
            // (U_ & fvc::reconstruct(fvc::snGrad(p_)* mesh_.magSf())) // Approach 2
            // (U_ & pGradLS_) // Approach 3
        )
    );

    volScalarField& pGradKE = tPGradKE.ref();
    pGradKE.write(runTime_.outputTime());

    // volume weighted average
    avgGradpKE_ = pGradKE.weightedAverage(mesh_.V());
    Info << "\tavg of pGrad KE: " << avgGradpKE_.value() << endl;

    return tPGradKE;
}


tmp<volScalarField> Foam::kineticEnergyAnalysis::getDissipationKE()
{
    // Calculate dissipation
    tmp<volScalarField> nuEff_
    (
        mesh_.lookupObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        ).nuEff()
    );

    tmp<volScalarField> tDissipKE
    (
        new volScalarField
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
        )
    );
    nuEff_.clear();

    volScalarField& dissipKE = tDissipKE.ref();
    dissipKE.write(runTime_.outputTime());

    // volume weighted average
    avgDisspKE_ = dissipKE.weightedAverage(mesh_.V());
    Info << "\tavg of dissipation KE: " << avgDisspKE_.value() << endl;

    return tDissipKE;
}

// public member functions ---------------------------------------------------//
void Foam::kineticEnergyAnalysis::analyzeKEBalance()
{

    tmp<volScalarField> tKEbalance
    (
        new volScalarField
        (
            IOobject
            (
                "KEbalance",
                U_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            // mesh_,
            // dimensionedScalar("", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0)
            (  getTemporalKE()
             + getConvectionKE()
             - getPGradKE()
             - getDissipationKE()
            )
        )
    );
    volScalarField& KEbalance = tKEbalance.ref();
    KEbalance.write(runTime_.outputTime());

    // avgKEBalance_ = (avgDdtKE_+avgConvKE_+avgGradpKE_+avgDisspKE_);
    avgKEBalance_ = KEbalance.weightedAverage(mesh_.V());
    Info << "\tavg of KE balance: " << avgKEBalance_.value() << endl;

    return;
}


// Error introduced in correcting U after poisson pressure correction of p & phi
// staggered gradient is used to correct face flux phi
// collocated gradient of pp is used to correct collocated velocity
// the term below is difference between to pressure gradient formulations
void Foam::kineticEnergyAnalysis::getPPGradDiffKE()
{
    const volScalarField& pp_ = mesh_.lookupObject<volScalarField>("pp");

    tmp<volScalarField> tppGradDiffKE
    (
        new volScalarField
        (
            IOobject
            (
                "ppGradDiffKE",
                pp_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            // mesh_,
            // dimensionedScalar("", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0)
            U_&
            ((fvc::grad(pp_)-fvc::reconstruct(fvc::snGrad(pp_)* mesh_.magSf())))
        )
    );

    volScalarField& ppGradDiffKE = tppGradDiffKE.ref();
    ppGradDiffKE.write(runTime_.outputTime());

    // volume weighted average
    avgPpGradKE_ = ppGradDiffKE.weightedAverage(mesh_.V());
    Info << "\tavg of ppGradDiff KE: " << avgPpGradKE_.value() << endl;

    return;
}


// the total contribution of the divegence error of collocated velocity
void Foam::kineticEnergyAnalysis::getPDivU()
{
    const volScalarField& pp_ = mesh_.lookupObject<volScalarField>("pp");

    tmp<volScalarField> tgradPpKE
    (
        new volScalarField
        (
            IOobject
            (
                "gradPpKE",
                pp_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            // pp_*fvc::div(U_)
            p_*fvc::div(U_)
        )
    );

    volScalarField& gradPpKE = tgradPpKE.ref();
    gradPpKE.write(runTime_.outputTime());

    // volume weighted average
    avgPDivU_ = gradPpKE.weightedAverage(mesh_.V());
    Info << "\tavg of gradPpKE KE: " << avgPDivU_.value() << endl;

    return;
}

void Foam::kineticEnergyAnalysis::getAddtionalKETerms()
{
    // pressure grad difference KE
    getPPGradDiffKE();

    // collocated velocity divergence error KE
    getPDivU();

    return;
}



void Foam::kineticEnergyAnalysis::setPropertiesOutput()
{
    if(Pstream::parRun() && !(Pstream::master()))
    {
        return;
    }

    // create output file
    fileName outputDir;
    autoPtr<fileName> outFilePath_;

    if(Pstream::parRun() && Pstream::master())
    {
        outputDir = runTime_.path()/"../postProcessing";

        outFilePath_.reset( new fileName
        (
            runTime_.path()/"../postProcessing"/
            ("kineticEnergyAnalysis.dat_"+runTime_.timeName())
        )
        );

    }
    else
    {
        outputDir = runTime_.path()/"postProcessing";

        outFilePath_.reset( new fileName
        (
            runTime_.path()/"postProcessing"/
            ("kineticEnergyAnalysis.dat_"+runTime_.timeName())
        )
        );
    }

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    KEPropertiesFile_.reset(new OFstream(*outFilePath_));

    // write header of the log file
    KEPropertiesFile_()
        << "time" << tab
        << "avgKEBalance(t)" << tab
        << "avgDdtKE(t)" << tab
        << "avgConvKE(t)" << tab
        << "avgGradpKE(t)" << tab
        << "avgDisspKE(t)" << tab
        << "avgPpGradKE(t)" << tab
        << "avgPDivU(t)" << endl;

    KEPropertiesFile_().precision(12); // set precision

    Info << "KE analysis: Global properties are written in\n"
         << KEPropertiesFile_().name() << endl;

    return;
}


void Foam::kineticEnergyAnalysis::writeAvgValues()
{
    if(Pstream::parRun() && !(Pstream::master()))
    {
        return;
    }

    Info << "Writing to KE-log file" << endl;

    KEPropertiesFile_()
        << runTime_.value() << tab
        << avgKEBalance_.value() << tab
        << avgDdtKE_.value() << tab
        << avgConvKE_.value() << tab
        << avgGradpKE_.value() << tab
        << avgDisspKE_.value() << tab
        << avgPpGradKE_.value() << tab
        << avgPDivU_.value() << endl;

    return;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticEnergyAnalysis::kineticEnergyAnalysis
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& p
)
:
    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    // Set the pointer to the velocity, flux and pressure-correction field
    U_(U),
    phi_(phi),
    p_(p),

    avgKEBalance_("avgKE",       dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0),
    avgDdtKE_    ("avgDdtKE",    dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0),
    avgConvKE_   ("avgConvKE",   dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0),
    avgGradpKE_  ("avgGradpKE",  dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0),
    avgDisspKE_  ("avgDisspKE",  dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0),
    avgPpGradKE_ ("avgPpGradKE", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0),
    avgPDivU_    ("avgPDivU", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0)
{
    // setPropertiesOutput();
    // analyzeKEBalance();
    // getPPGradDiffKE();
    // writeAvgValues();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticEnergyAnalysis::~kineticEnergyAnalysis()
{
    Info << "kineticEnergyAnalysis Destructor" << endl;
}


// ************************************************************************* //
