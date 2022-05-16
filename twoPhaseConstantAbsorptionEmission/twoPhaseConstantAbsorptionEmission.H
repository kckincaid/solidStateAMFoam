/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::radiation::twoPhaseConstantAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description
    Constant radiation absorption and emission coefficients for two phase
	fluid

SourceFiles
    twoPhaseConstantAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_twoPhaseConstantAbsorptionEmission_H
#define radiation_twoPhaseConstantAbsorptionEmission_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class twoPhaseConstantAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseConstantAbsorptionEmission
:
    public absorptionEmissionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar a1_;
		dimensionedScalar a2_;

        //- Emission coefficient / [1/m]
        dimensionedScalar e1_;
		dimensionedScalar e2_;

        //- Emission contribution / [kg/(m s^3)]
        dimensionedScalar E1_;
		dimensionedScalar E2_;

		//- Phase fraction field to read in
		word alphaField_;

		// Phase fraction used to calculate average values
		//volScalarField alpha_;


public:

    //- Runtime type information
    TypeName("twoPhaseConstantAbsorptionEmission");


    // Constructors

        //- Construct from components
        twoPhaseConstantAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~twoPhaseConstantAbsorptionEmission();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
