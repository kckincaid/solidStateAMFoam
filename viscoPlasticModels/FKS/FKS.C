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

#include "FKS.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(FKS, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        FKS,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

// Calculate the strain hardening contribution
//Foam::tmp<Foam::volScalarField>
//Foam::viscosityModels::FKS::H() const
//{
//	// Import strain rate and plastic strain fields
//	const volScalarField& edot_ = U_.mesh().lookupObject<volScalarField>("epsilonDotEq");
//	const volScalarField& ep_ = U_.mesh().lookupObject<volScalarField>("epsilonEq");

//    return
//	(
//		a1_*(edot_/edot_) + a2_*atan(a3_*ep_)
//	);
//}
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::FKS::H() const
{
	// Import strain rate field
	const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");

    return
	(
		// Return a scalar field equivalent to the strain hardening contribution at full
		// saturation. Note that atan(x) --> pi/2 when x >> 1
		(a1_ + a2_*1.571)*temp_/temp_
	);
}

// Calculate the strain rate stiffening contribution
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::FKS::Lambda(const volScalarField Th_) const
{
	// Import strain rate field
	const volScalarField& edot_ = U_.mesh().lookupObject<volScalarField>("epsilonDotEq");

	// Create smol temp value
	dimensionedScalar Ts_("Ts", dimless, SMALL);

    return
    (
		1.0 + (b1_*pow((Th_ + Ts_), b2_))*(b3_*Foam::log(edot_/e0_ + VSMALL))
    );
}

// Calculate the thermal softening contribution
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::FKS::Theta(const volScalarField Th_) const
{
    return
    (
		1.0 - 1.0/pow(1.0 + exp(-c1_*Th_), 1/c2_)
    );
}

// Calculate the homologous temperature
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::FKS::Tstar(const volScalarField T_) const
{
    return
    (
		(T_ - Ta_)/(Tm_ - Ta_)
    );
}

// Calculate the flow stress
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::FKS::sigmaf() const
{
	// Import temperature field
	const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");

	// Convert to homologous temperature
	volScalarField Th_ = Tstar(temp_);

	// Return flow stress
    return
    (
		H()*Lambda(Th_)*Theta(Th_)
    );
}


// Calculate effective kinematic viscosity
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::FKS::calcNu() const
{
	// Import effective strain rate field
	const volScalarField& edot_ = U_.mesh().lookupObject<volScalarField>("epsilonDotEq");

	// Look up limiters or set to default value if not present
	dimensionedScalar nuMin_("nuMin", dimViscosity, FKSCoeffs_.lookupOrDefault("nuMin", ROOTVSMALL));
	dimensionedScalar nuMax_("nuMax", dimViscosity, FKSCoeffs_.lookupOrDefault("nuMax", ROOTVGREAT));
	dimensionedScalar strainRateEffMin_("strainRateEffMin", dimless/dimTime, FKSCoeffs_.lookupOrDefault("strainRateEffMin", ROOTVSMALL));

	// Calculate effective viscosity field
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
			sigmaf()/(3.0*rho_*max(edot_,strainRateEffMin_))
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::FKS::FKS
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    FKSCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),

    a1_("a1", dimPressure, FKSCoeffs_),
	a2_("a2", dimPressure, FKSCoeffs_),
	//a3_("a3", dimless, FKSCoeffs_),

	b1_("b1", dimless, FKSCoeffs_),
	b2_("b2", dimless, FKSCoeffs_),
	b3_("b3", dimless, FKSCoeffs_),
	e0_("e0", dimless/dimTime, FKSCoeffs_),

	c1_("c1", dimless, FKSCoeffs_),
	c2_("c2", dimless, FKSCoeffs_),
	Tm_("Tm", dimTemperature, FKSCoeffs_),
	Ta_("Ta", dimTemperature, FKSCoeffs_),

	rho_("rho", dimDensity, FKSCoeffs_),

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
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::FKS::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    FKSCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    FKSCoeffs_.lookup("a1") >> a1_;
    FKSCoeffs_.lookup("a2") >> a2_;
    //FKSCoeffs_.lookup("a3") >> a3_;

    FKSCoeffs_.lookup("b1") >> b1_;
    FKSCoeffs_.lookup("b2") >> b2_;
    FKSCoeffs_.lookup("b3") >> b3_;
    FKSCoeffs_.lookup("e0") >> e0_;

    FKSCoeffs_.lookup("c1") >> c1_;
    FKSCoeffs_.lookup("c2") >> c2_;
    FKSCoeffs_.lookup("Tm") >> Tm_;
    FKSCoeffs_.lookup("Ta") >> Ta_;

	FKSCoeffs_.lookup("rho") >> rho_;

    return true;
}

// ************************************************************************* //
