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
    Foam::viscosityModels::FKS

Description
    Modified Fraser-Kiss-St. Georges viscoplastic flow model. The original
	model included a strain-hardening function which is set to unity here
	due to the limitations of strain tracking with FVM. Requires ten (10)
	material-specific constants which may be available in the literature.

SourceFiles
    FKS.C

\*---------------------------------------------------------------------------*/

#ifndef FKS_H
#define FKS_H

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
                           Class FKS Declaration
\*---------------------------------------------------------------------------*/

class FKS
:
    public viscosityModel
{
    // Private data

        dictionary FKSCoeffs_;

		// Strain hardening constants
		dimensionedScalar a1_;
		dimensionedScalar a2_;
		//dimensionedScalar a3_;

		// Strain rate stiffening constants
        dimensionedScalar b1_;
		dimensionedScalar b2_;
		dimensionedScalar b3_;
		dimensionedScalar e0_;

		// Thermal softening constants
		dimensionedScalar c1_;
		dimensionedScalar c2_;
		dimensionedScalar Tm_;
		dimensionedScalar Ta_;

		// Material density
		dimensionedScalar rho_;

protected:

    // Protected data

        volScalarField nu_;
		volScalarField sigmaf_;
		volScalarField sigmay_;

    // Protected Member Functions

		//- Calculate and return the strain hardening function
		//tmp<volScalarField> H() const;
		tmp<volScalarField> H() const;

		//- Calculate and return the strain rate stiffening function
		tmp<volScalarField> Lambda(const volScalarField Th_) const;

		//- Calculate and return the thermal softening function
		tmp<volScalarField> Theta(const volScalarField Th_) const;

		//- Calculate and return the homologous temperature
		tmp<volScalarField> Tstar(const volScalarField T_) const;

		//- Calculate and return the flow stress
		tmp<volScalarField> sigmaf() const;

		//- Calculate and return the yield stress (at 0 plastic strain, 0.1 strain rate)
		tmp<volScalarField> sigmay() const;

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("FKS");


    // Constructors

        //- Construct from components
        FKS
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~FKS()
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
			sigmaf_ = sigmaf();
			sigmay_ = sigmay();         
			nu_ = calcNu();
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
