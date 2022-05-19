/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "slipStick.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace stickModels
{
    defineTypeNameAndDebug(slipStick, 0);
    addToRunTimeSelectionTable
    (
        stickModel,
        slipStick,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

// Calculate viscous heating at cell centers
Foam::tmp<Foam::volScalarField>
Foam::stickModels::slipStick::calcQvisc() const
{
	// Import viscosity and effective strain rate fields
	const volScalarField& muEff_ = U_.mesh().lookupObject<volScalarField>("muEff");
	const volScalarField& eDotEq_ = U_.mesh().lookupObject<volScalarField>("epsilonDotEq");

	// Calculate and return viscous heating
	//qvisc_ = 3.0*phiVisc_*muEff_*eDotEq_*eDotEq_;

    //return qvisc_;
	return(3.0*phiVisc_*muEff_*eDotEq_*eDotEq_);
}

// Calculate frictional heating on patch
Foam::tmp<Foam::surfaceScalarField>
Foam::stickModels::slipStick::calcQfric() const
{
	// Create face-interpolated phase fraction and pressure
	surfaceScalarField alphaf_
    (
        fvc::interpolate(min(max(alpha_, scalar(0)), scalar(1)))
    );

	surfaceScalarField pf_
    (
        fvc::interpolate(p_)
    );

	// Create unity value scalar field
	surfaceScalarField one_ = alphaf_/(alphaf_+VSMALL);

	//qfric_ = alphaf_*((one_-delta())*etaf_*tau() + delta()*muf()*pf_)
	  	   //* (mag(omega_)*r()); // add in Vt sin theta contribution

	// Calculate and return frictional heating field
    //return qfric_;
	return
	(
		alphaf_*((one_-delta())*etaf_*tau() + delta()*muf()*pf_)
	  * (mag(omega_)*mag(r())) // add in Vt sin theta contribution
	);
}

// Calculate distance from center of patch
Foam::tmp<Foam::surfaceVectorField>
Foam::stickModels::slipStick::r() const
{
	// Get patch info
	const scalar patchNum_ = U_.mesh().boundaryMesh().findPatchID(patch_);
	const polyPatch& currPatch_ = U_.mesh().boundaryMesh()[patchNum_];

	// Initialize center vector and area
	vector patchCenter_ = Zero;
	scalar patchArea_ = VSMALL;

	// Create pointers to face surface area and centers (used by both loops)
	const surfaceScalarField& magSf_ = U_.mesh().magSf();
	const surfaceVectorField& Cf_ = U_.mesh().Cf();

	// Calculate area-weighted patch center
	forAll(currPatch_, facei_)
	{
		patchCenter_ += Cf_.boundaryField()[patchNum_][facei_]*magSf_.boundaryField()[patchNum_][facei_];
		patchArea_ += magSf_.boundaryField()[patchNum_][facei_];
	}

	// Sum over all processors
	reduce(patchCenter_, sumOp<vector>());
	reduce(patchArea_, sumOp<scalar>());

	// Calculate patch centroid
	patchCenter_ = patchCenter_ / patchArea_;

	// Create surfaceScalarField of patch center
	surfaceVectorField c_(patchCenter_*magSf_/magSf_);

	// Return info on patch center location
	Info << "Center for patch " << patch_ << " found to be at " << patchCenter_ << endl;

    return(Cf_ - c_);
}

// Calculate slip fraction
Foam::tmp<Foam::surfaceScalarField>
Foam::stickModels::slipStick::delta() const
{
	// Calculate patch velocity
	surfaceVectorField Up_
	(
		omega_^r()
	);

	// Calculate slip fraction bsed on actual velocity
	const surfaceVectorField Uf_ = fvc::interpolate(U_);

	surfaceScalarField delta_ = mag(Up_ - Uf_);

    return delta_;
}

// Calculate friction coefficient
Foam::tmp<Foam::surfaceScalarField>
Foam::stickModels::slipStick::muf() const
{
	return
	(
		mu0_*exp(-lambda_*delta()*mag(omega_)*mag(r()))
	);
}

// Calculate temp-dependent max shear stress
Foam::tmp<Foam::surfaceScalarField>
Foam::stickModels::slipStick::tau() const
{
	// Import and limit temperature
	const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");

	const surfaceScalarField tempLim_
	(
		fvc::interpolate(min(max(temp_, Tmin_), Tmax_))
	);

	// Calculate and return max shear stress
	return
	(
		TC0_*tempLim_/tempLim_ + TC1_*tempLim_ + TC2_*pow(tempLim_, 2)
	  + TC3_*pow(tempLim_, 3) + TC4_*pow(tempLim_, 4)
	);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stickModels::slipStick::slipStick
(
    const word& name,
    const dictionary& stickProperties,
    const volVectorField& U,
    const volScalarField& alpha,
	const volScalarField& p
)
:
    stickModel(name, stickProperties, U, alpha, p),
    slipStickCoeffs_(stickProperties.optionalSubDict(typeName + "Coeffs")),

	phiVisc_("phiVisc", dimless, slipStickCoeffs_),

	patch_(stickProperties.lookup("fricPatch")),
    mu0_("mu0", dimless, slipStickCoeffs_),
    lambda_("lambda", dimTime/dimLength, slipStickCoeffs_),
    omega_("omega", dimless/dimTime, slipStickCoeffs_),
    etaf_("etaf", dimless, slipStickCoeffs_),

    TC0_("TC0", dimPressure, slipStickCoeffs_),
    TC1_("TC1", dimPressure/dimTemperature, slipStickCoeffs_),
    TC2_("TC2", dimPressure/pow(dimTemperature,2), slipStickCoeffs_),
    TC3_("TC3", dimPressure/pow(dimTemperature,3), slipStickCoeffs_),
    TC4_("TC4", dimPressure/pow(dimTemperature,4), slipStickCoeffs_),
	Tmax_("Tmax", dimTemperature, slipStickCoeffs_),
	Tmin_("Tmin", dimTemperature, slipStickCoeffs_),

    qvisc_
    (
        IOobject
        (
            "qvisc",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    calcQvisc()
    ),

    qfric_
    (
        IOobject
        (
            "qfric",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    calcQfric()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::stickModels::slipStick::read
(
    const dictionary& stickProperties
)
{
    stickModel::read(stickProperties);

    slipStickCoeffs_ = stickProperties.optionalSubDict(typeName + "Coeffs");

	slipStickCoeffs_.lookup("phiVisc") >> phiVisc_;

	slipStickCoeffs_.lookup("fricPatch") >> patch_;
    slipStickCoeffs_.lookup("mu0") >> mu0_;
    slipStickCoeffs_.lookup("lambda") >> lambda_;
    slipStickCoeffs_.lookup("omega") >> omega_;
    slipStickCoeffs_.lookup("etaf") >> etaf_;

	slipStickCoeffs_.lookup("TC0") >> TC0_;
	slipStickCoeffs_.lookup("TC1") >> TC1_;
	slipStickCoeffs_.lookup("TC2") >> TC2_;
	slipStickCoeffs_.lookup("TC3") >> TC3_;
	slipStickCoeffs_.lookup("TC4") >> TC4_;
	slipStickCoeffs_.lookup("Tmax") >> Tmax_;
	slipStickCoeffs_.lookup("Tmin") >> Tmin_;

    return true;
}


// ************************************************************************* //
