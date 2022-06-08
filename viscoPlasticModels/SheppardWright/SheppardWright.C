/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "SheppardWright.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(SheppardWright, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        SheppardWright,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

// Calculate Zener-Hollomon parameter
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::SheppardWright::Z() const
{
	// Import temperature field
	const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");

	// Import effective strain rate field
	const volScalarField& epsilonDotEq_ = U_.mesh().lookupObject<volScalarField>("epsilonDotEq");

	// Return Z
    return
    (
		//sqrt(2.0/3.0)*strainRate()*exp(Q_/(R_*temp_))
		epsilonDotEq_*exp(Q_/(R_*temp_))
    );
}

// Calculate flow stress
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::SheppardWright::sigmaf() const
{
    return
    (
		sigmaR_*asinh(pow(Z()/A_, 1/n_))
    );
}

// Calculate the yield stress (at a strain rate of 0.1 1/s)
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::SheppardWright::sigmay() const
{
	// Import temperature field
	const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");

	// Create dimensioned scalar with value of 0.1 1/s
	dimensionedScalar e01("e01", dimless/dimTime, 0.1);

	// Return yield stress at given temp and strain rate of 0.1
    return
    (
		sigmaR_*asinh(pow(e01*exp(Q_/(R_*temp_))/A_, 1/n_))
    );
}


// Calculate effective kinematic viscosity
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::SheppardWright::calcNu() const
{
	// Import effective strain rate field
	const volScalarField& epsilonDotEq_ = U_.mesh().lookupObject<volScalarField>("epsilonDotEq");

	// Look up limiters or set to default value if not present
	dimensionedScalar nuMin_("nuMin", dimViscosity, SheppardWrightCoeffs_.lookupOrDefault("nuMin", ROOTVSMALL));
	dimensionedScalar nuMax_("nuMax", dimViscosity, SheppardWrightCoeffs_.lookupOrDefault("nuMax", ROOTVGREAT));
	dimensionedScalar eDotMin_("eDotMin", dimless/dimTime, SheppardWrightCoeffs_.lookupOrDefault("eDotMin", ROOTVSMALL));

	// Calculate effective viscosity field
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
			sigmaf()/(3.0*rho_*max(epsilonDotEq_, eDotMin_))
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::SheppardWright::SheppardWright
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    SheppardWrightCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    sigmaR_("sigmaR", dimPressure, SheppardWrightCoeffs_),
	Q_("Q", dimEnergy/dimMoles, SheppardWrightCoeffs_),
	R_("R", dimEnergy/dimMoles/dimTemperature, SheppardWrightCoeffs_),
	A_("A", dimless/dimTime, SheppardWrightCoeffs_),
    n_("n", dimless, SheppardWrightCoeffs_),
	rho_("rho", dimDensity, SheppardWrightCoeffs_),

    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    ),

    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	sigmaf()
    ),

    sigmay_
    (
        IOobject
        (
            "sigmay",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	sigmay()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::SheppardWright::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    SheppardWrightCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    SheppardWrightCoeffs_.lookup("sigmaR") >> sigmaR_;
    SheppardWrightCoeffs_.lookup("Q") >> Q_;
	SheppardWrightCoeffs_.lookup("R") >> R_;
	SheppardWrightCoeffs_.lookup("A") >> A_;
	SheppardWrightCoeffs_.lookup("n") >> n_;

    return true;
}

// ************************************************************************* //
