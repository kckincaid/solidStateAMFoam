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

    // Return viscous heating
	return(3.0*phiVisc_*muEff_*pow(eDotEq_, 2));
}

// Calculate frictional heating on patch
Foam::tmp<Foam::volScalarField>
Foam::stickModels::slipStick::calcQfric() const
{
	// Create vector field of distance from patch center
	volVectorField rad_ = r();

	// Create a volScalarField with a value of 1 for cells adjacent to friction
	// patch and a value of 0 for all other cells
	volScalarField fricCell_ = alpha_*0.0;

	// Get patch info
	scalar patchNum_ = U_.mesh().boundaryMesh().findPatchID(patch_);
	const polyPatch& currPatch_ = U_.mesh().boundaryMesh()[patchNum_];

	// Loop over faces on patch
	forAll(currPatch_, facei_)
	{
		label faceCelli_ = currPatch_.faceCells()[facei_];

		// Set fricCell to 1.0 in cells adjacent to patch
		fricCell_[faceCelli_] = 1.0;
	}

	dimensionedScalar ASMALL("VVSMALL", dimVelocity*dimLength, SMALL);

	// Get updated alpha field (one from constructor doesn't update)
	const volScalarField& alpha1_ = U_.mesh().lookupObject<volScalarField>(fricPhaseName_);

	// Return frictional heating
	return
	(
		alpha1_*((scalar(1.0)-delta(rad_))*etaf_*tau() + delta(rad_)*muf(rad_)*p_)
	  * (mag(omega_)*mag(rad_) - mag(Ut_)*Foam::sqrt(1.0 - sqr((rad_ & Ut_)/(mag(rad_)*mag(Ut_)+ASMALL))))
	);
}

// Calculate distance from center of patch as a vector field
Foam::tmp<Foam::volVectorField>
Foam::stickModels::slipStick::r() const
{
	// Get info for patch
	scalar patchNum_ = U_.mesh().boundaryMesh().findPatchID(patch_);
	const polyPatch& currPatch_ = U_.mesh().boundaryMesh()[patchNum_];

	// Initialize center vector and area
	vector patchCenter_ = Zero;
	scalar patchArea_ = VSMALL;

	// Create pointers to face surface area and centers
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

	// Create a dimensioned vector from patchCenter vector
	dimensionedVector C_("C", dimLength, patchCenter_);

	// Return distance from patch centroid
	return(U_.mesh().C() - C_);	
}

// Calculate slip fraction as a volScalarField
Foam::tmp<Foam::volScalarField>
Foam::stickModels::slipStick::delta(const volVectorField r_) const
{
	//return(mag(U_)/mag(r_^omega_));
	return
	(
		scalar(1.0) - exp(-mag(omega_)*mag(r_)/(delta0_*omega0_*Rs_))
	);	
}

// Calculate friction coefficient as a volScalarField
Foam::tmp<Foam::volScalarField>
Foam::stickModels::slipStick::muf(const volVectorField r_) const
{
	return(mu0_*exp(-lambda_*delta(r_)*mag(omega_)*mag(r_)));	
}

// Calculate maximum shear stress as a volScalarField
Foam::tmp<Foam::volScalarField>
Foam::stickModels::slipStick::tau() const
{
	// Create pointer to temp field and create limited temp
	const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");

	const volScalarField TLim_
	(
		min(max(temp_, Tmin_), Tmax_)
	);

	return
	(
		(TC0_*TLim_/TLim_ + TC1_*TLim_ + TC2_*pow(TLim_, 2) 
	  + TC3_*pow(TLim_, 3) + TC4_*pow(TLim_, 4))*1.0E+6
	);	
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

	fricPhaseName_(stickProperties.lookup("fricPhaseName")),

	patch_(stickProperties.lookup("fricPatch")),
    mu0_("mu0", dimless, slipStickCoeffs_),
    lambda_("lambda", dimTime/dimLength, slipStickCoeffs_),
    omega_("omega", dimless/dimTime, slipStickCoeffs_),
    etaf_("etaf", dimless, slipStickCoeffs_),

	delta0_("delta0", dimless, slipStickCoeffs_),
	omega0_("omega0", dimless/dimTime, slipStickCoeffs_),
	Rs_("Rs", dimLength, slipStickCoeffs_),

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
            U_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    U_.mesh()
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

	slipStickCoeffs_.lookup("fricPhaseName") >> fricPhaseName_;
	slipStickCoeffs_.lookup("fricPatch") >> patch_;
    slipStickCoeffs_.lookup("mu0") >> mu0_;
    slipStickCoeffs_.lookup("lambda") >> lambda_;
    slipStickCoeffs_.lookup("omega") >> omega_;
    slipStickCoeffs_.lookup("etaf") >> etaf_;

	slipStickCoeffs_.lookup("delta0") >> delta0_;
	slipStickCoeffs_.lookup("omega0") >> omega0_;
	slipStickCoeffs_.lookup("Rs") >> Rs_;

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
