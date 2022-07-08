/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd
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
    Foam::stickModels::exponentialSlipStick

Description
    A slip-stick model in which the slip fraction and friction coefficient
	vary as a function of radius, or distance from the tool center. Requires
	a number of user-defined inputs which control the degree of slip and 
	magnitude of friction coefficient.

SourceFiles
    exponentialSlipStick.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialSlipStick_H
#define exponentialSlipStick_H

#include "stickModel.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stickModels
{


/*---------------------------------------------------------------------------*\
                           Class exponentialSlipStick Declaration
\*---------------------------------------------------------------------------*/

class exponentialSlipStick
:
    public stickModel
{
    // Private data

        dictionary exponentialSlipStickCoeffs_;

		// Names of fields and patches
		word alphaName_;
		word pName_;
		word patch_;

		// Coefficients for calculations
		dimensionedScalar betaTQ_;	// Taylor-Quinney coefficient
        dimensionedVector omega_;	// Angular velocity vector
		dimensionedScalar omega0_;	// Slip model coefficient #1
		dimensionedScalar delta0_;	// Slip model coefficient #2
		dimensionedScalar Rs_;		// Tool shoulder radius
        dimensionedScalar mu0_;		// Friction model coefficient #1
        dimensionedScalar lambda_;	// Friction model coefficient #2

protected:

    // Protected data

        volScalarField qvisc_;
		volScalarField qfric_;


    // Protected Member Functions

        //- Calculate and return the heating terms
        tmp<volScalarField> calcQvisc() const;
		tmp<volScalarField> calcQfric() const;

		//- Other functions used in calculations
		tmp<volVectorField> r()	const;
		tmp<volScalarField> delta(const volVectorField r_) const;
		tmp<volScalarField> muf(const volVectorField r_) const;


public:

    //- Runtime type information
    TypeName("exponentialSlipStick");


    // Constructors

        //- Construct from components
        exponentialSlipStick
        (
            const word& name,
            const dictionary& stickProperties,
            const volVectorField& U,
			const volVectorField& Ut
        );


    //- Destructor
    virtual ~exponentialSlipStick()
    {}


    // Member Functions

        //- Return the cell-centered viscous heating
        tmp<volScalarField> qvisc() const
        {
            return qvisc_;
        }

        //- Return the cell-centered viscous heating
        tmp<volScalarField> qfric() const
        {
            return qfric_;
        }

        //- Correct the heating terms
        virtual void correct()
        {
            qvisc_ = calcQvisc();
			qfric_ = calcQfric();
			qfric_.correctBoundaryConditions();
        }

        //- Read frictionProperties dictionary
        virtual bool read(const dictionary& stickProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stickModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
