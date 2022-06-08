/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::viscosityModels::SheppardWright

Description
    Sheppard-Wright viscoplastic model. Calculates viscosity as a function
	of strain rate and temperature. Requires six material-specific model
	constants to be specified by the user, which are widely available in the
	literature for many metals.

SourceFiles
    SheppardWright.C

\*---------------------------------------------------------------------------*/

#ifndef SheppardWright_H
#define SheppardWright_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class SheppardWright Declaration
\*---------------------------------------------------------------------------*/

class SheppardWright
:
    public viscosityModel
{
    // Private data

        dictionary SheppardWrightCoeffs_;

        dimensionedScalar sigmaR_;
        dimensionedScalar Q_;
		dimensionedScalar R_;
		dimensionedScalar A_;
		dimensionedScalar n_;
		dimensionedScalar rho_;

protected:

    // Protected data

        volScalarField nu_;
		volScalarField sigmaf_;	// Flow stress field
		volScalarField sigmay_; // Yield stress field for qfric calc

    // Protected Member Functions

		//- Calculate and return the Zener-Hollomon parameter
		tmp<volScalarField> Z() const;

		//- Calculate and return the flow stress
		tmp<volScalarField> sigmaf() const;

		//- Calculate and return the yield stress
		tmp<volScalarField> sigmay() const;

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("SheppardWright");


    // Constructors

        //- Construct from components
        SheppardWright
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~SheppardWright()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();
			sigmaf_ = sigmaf();
			sigmay_ = sigmay();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
