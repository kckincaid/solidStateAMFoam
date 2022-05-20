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
    Foam::stickModels::slipStick

Description
    A basic slip-stick model which works with an arbitary slip boundary condition
	(e.g. constant, exponential, etc.) and a constant or exponential friction
	coefficient.

SourceFiles
    slipStick.C

\*---------------------------------------------------------------------------*/

#ifndef slipStick_H
#define slipStick_H

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
                           Class slipStick Declaration
\*---------------------------------------------------------------------------*/

class slipStick
:
    public stickModel
{
    // Private data

        dictionary slipStickCoeffs_;

		// Viscous heating premultiplier
		dimensionedScalar phiVisc_;

		// Friction quantities
		word patch_;
        dimensionedScalar mu0_;
        dimensionedScalar lambda_;
        dimensionedVector omega_;
        dimensionedScalar etaf_;

		// Temp-dependent max shear stress quantities
		dimensionedScalar TC0_;
		dimensionedScalar TC1_;
		dimensionedScalar TC2_;
		dimensionedScalar TC3_;
		dimensionedScalar TC4_;
		dimensionedScalar Tmax_;
		dimensionedScalar Tmin_;


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
		tmp<volScalarField> tau() const;


public:

    //- Runtime type information
    TypeName("slipStick");


    // Constructors

        //- Construct from components
        slipStick
        (
            const word& name,
            const dictionary& stickProperties,
            const volVectorField& U,
			const volVectorField& Ut,
            const volScalarField& alpha,
			const volScalarField& p
        );


    //- Destructor
    virtual ~slipStick()
    {}


    // Member Functions

        //- Return the cell-centered viscous heating
        virtual tmp<volScalarField> qvisc() const
        {
            return qvisc_;
        }

        //- Return the cell-centered viscous heating
        virtual tmp<volScalarField> qfric() const
        {
            return qfric_;
        }

        //- Correct the heating terms
        virtual void correct()
        {
            qvisc_ = calcQvisc();
			qfric_ = calcQfric();
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
