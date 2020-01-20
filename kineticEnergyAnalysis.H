/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::kineticEnergyAnalysis

Description
    Term-by-term kinetic energy analysis object

SourceFiles
    kineticEnergyAnalysis.C

\*---------------------------------------------------------------------------*/

#ifndef kineticEnergyAnalysis_H
#define kineticEnergyAnalysis_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class kineticEnergyAnalysis Declaration
\*---------------------------------------------------------------------------*/

class kineticEnergyAnalysis
{
    // Private data

        //- Runtime pointer
        const Time& runTime_;;

        //- Mesh pointer
        const fvMesh& mesh_;

        // Pointer to solver fields
        //velocity field pointer
        const volVectorField& U_;

        // velocity flux field pointer
        const surfaceScalarField& phi_;

        // pressure field pointer
        const volScalarField& p_;

    // Private member function
        tmp<volScalarField> getTemporalKE();
        tmp<volScalarField> getConvectionKE();
        tmp<volScalarField> getPGradKE();
        tmp<volScalarField> getDissipationKE();


public:

    //- Declare name of the class and its debug switch
    ClassName("kineticEnergyAnalysis");


    // Constructors
    kineticEnergyAnalysis
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& p
    );


    // Destructor
    virtual ~kineticEnergyAnalysis();


    // Public Member functions
        void analyzeKEBalance();

        void getPPGradDiffKE();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
