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
	// Create pointers to viscosity and effective strain rate fields
	const volScalarField& muEff_ = U_.mesh().lookupObject<volScalarField>("muEff");
	const volScalarField& eDotEq_ = U_.mesh().lookupObject<volScalarField>("epsilonDotEq");

    // Return qvisc_;
	return(3.0*phiVisc_*muEff_*pow(eDotEq_, 2));
}

// Calculate frictional heating on patch
Foam::tmp<Foam::volScalarField>
Foam::stickModels::slipStick::calcQfric() const
{
	// Import and limit temperature field for max shear stress
	const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");

	const volScalarField TLim_
	(
		min(max(temp_, Tmin_), Tmax_)
	);

	// Get patch info for qfric
	scalar patchNum_ = U_.mesh().boundaryMesh().findPatchID(patch_);
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
		// Add area-weighted face center and face area to totals
		patchCenter_ += Cf_.boundaryField()[patchNum_][facei_]*magSf_.boundaryField()[patchNum_][facei_];
		patchArea_ += magSf_.boundaryField()[patchNum_][facei_];
	}

	// Perform calculation for all processors in decomposed case
	reduce(patchCenter_, sumOp<vector>());
	reduce(patchArea_, sumOp<scalar>());

	// Calculate patch centroid
	patchCenter_ = patchCenter_ / patchArea_;

	Info << "Center for patch " << patch_ << " found to be at " << patchCenter_ << endl;

	// Loop over all faces on patch		
	forAll(currPatch_, facei_)
	{
		// Get labels for cells which own patch faces
		label faceCelli_ = currPatch_.faceCells()[facei_];

		// Get distance of cell center from center of rotation
		dimensionedScalar dx_("dx", dimLength, Cf_.boundaryField()[patchNum_][facei_].component(0) - patchCenter_.component(0));
		dimensionedScalar dy_("dy", dimLength, Cf_.boundaryField()[patchNum_][facei_].component(1) - patchCenter_.component(1));
		dimensionedScalar dz_("dz", dimLength, Cf_.boundaryField()[patchNum_][facei_].component(2) - patchCenter_.component(2));

		dimensionedScalar rad_ = Foam::sqrt(pow(dx_,2) + pow(dy_,2) + pow(dz_,2));

		// Calculate angle of point w.r.t origin
		scalar theta_ = Foam::atan(Cf_.boundaryField()[patchNum_][facei_].component(1)
					  / Cf_.boundaryField()[patchNum_][facei_].component(0));

		// Calculate local slip fraction based on patch velocity and actual velocity
		dimensionedScalar Up_ = rad_*mag(omega_);
		dimensionedScalar delta_("delta", dimless, mag(U_[faceCelli_])/Up_.value());

		// Calculate local friction coefficient
		dimensionedScalar muf_ = mu0_*exp(-lambda_*delta_*mag(omega_)*rad_);

		// Calculate local maximum shear stress
		dimensionedScalar tau_ = TC0_;

		tau_ = TC0_ + TC1_*TLim_[faceCelli_] + TC2_*pow(TLim_[faceCelli_], 2) 
			 + TC3_*pow(TLim_[faceCelli_], 3) + TC4_*pow(TLim_[faceCelli_], 4);

		// Calculate local frictional heating
		qfric_[faceCelli_] = alpha_[faceCelli_]*((scalar(1) - delta_.value())*etaf_.value()*tau_.value() 
						   + delta_.value()*muf_.value()*p_[faceCelli_])
						   * (mag(omega_).value()*mag(Ut_[faceCelli_])*Foam::sin(theta_));
	}

	return(qfric_);	
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stickModels::slipStick::slipStick
(
    const word& name,
    const dictionary& stickProperties,
    const volVectorField& U,
    const volVectorField& Ut,
	const volScalarField& alpha,
	const volScalarField& p
)
:
    stickModel(name, stickProperties, U, Ut, alpha, p),
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
