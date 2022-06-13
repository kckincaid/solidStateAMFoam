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

#include "constantSlipStick.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace stickModels
{
    defineTypeNameAndDebug(constantSlipStick, 0);
    addToRunTimeSelectionTable
    (
        stickModel,
        constantSlipStick,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

// Calculate viscous heating at cell centers
Foam::tmp<Foam::volScalarField>
Foam::stickModels::constantSlipStick::calcQvisc() const
{
	// Create pointers to viscosity and effective strain rate fields
	const volScalarField& muEff_ = U_.mesh().lookupObject<volScalarField>("muEff");
	const volScalarField& eDotEq_ = U_.mesh().lookupObject<volScalarField>("epsilonDotEq");

    // Return viscous heating
	return(3.0*betaTQ_*muEff_*pow(eDotEq_, 2));
}

// Calculate frictional heating on patch
Foam::tmp<Foam::volScalarField>
Foam::stickModels::constantSlipStick::calcQfric() const
{
	// Create vector field of distance from patch center
	volVectorField rad_ = r();

	// Get updated alpha field (one from constructor doesn't update)
	const volScalarField& alpha_ = U_.mesh().lookupObject<volScalarField>(alphaName_);

	// Get yield stress field created by viscoplastic model
	const volScalarField& sigmay_ = U_.mesh().lookupObject<volScalarField>("sigmay");

	// Get updated pressure field
	const volScalarField& p_ = U_.mesh().lookupObject<volScalarField>(pName_);

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

	// Create a small field to prevent FPE issues in qfric calc
	//dimensionedScalar ASMALL("VVSMALL", dimVelocity*dimLength, SMALL);

	// Return frictional heating
	return
	(
		fricCell_*alpha_*((scalar(1.0)-delta_)*betaTQ_*sigmay_/Foam::sqrt(3.0) + delta_*muf_*p_)*(mag(omega_)*mag(rad_))
	);
//	return
//	(
//		alpha1_*((scalar(1.0)-delta(rad_))*etaf_*tau() + delta(rad_)*muf(rad_)*p_)
//	  * (mag(omega_)*mag(rad_) - mag(Ut_)*Foam::sqrt(1.0 - sqr((rad_ & Ut_)/(mag(rad_)*mag(Ut_)+ASMALL))))
//	);
}

// Calculate distance from center of patch as a vector field
Foam::tmp<Foam::volVectorField>
Foam::stickModels::constantSlipStick::r() const
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stickModels::constantSlipStick::constantSlipStick
(
    const word& name,
    const dictionary& stickProperties,
    const volVectorField& U,
    const volVectorField& Ut
)
:
    stickModel(name, stickProperties, U, Ut),
    constantSlipStickCoeffs_(stickProperties.optionalSubDict(typeName + "Coeffs")),

	alphaName_(stickProperties.lookup("alphaName")),
	pName_(stickProperties.lookup("pName")),
	patch_(stickProperties.lookup("fricPatch")),

	betaTQ_("betaTQ", dimless, constantSlipStickCoeffs_),
	delta_("delta", dimless, constantSlipStickCoeffs_),
    omega_("omega", dimless/dimTime, constantSlipStickCoeffs_),
    muf_("muf", dimless, constantSlipStickCoeffs_),

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

bool Foam::stickModels::constantSlipStick::read
(
    const dictionary& stickProperties
)
{
    stickModel::read(stickProperties);

    constantSlipStickCoeffs_ = stickProperties.optionalSubDict(typeName + "Coeffs");

	constantSlipStickCoeffs_.lookup("alphaName") >> alphaName_;
	constantSlipStickCoeffs_.lookup("pName") >> pName_;
	constantSlipStickCoeffs_.lookup("fricPatch") >> patch_;

	constantSlipStickCoeffs_.lookup("betaTQ") >> betaTQ_;
	constantSlipStickCoeffs_.lookup("delta") >> delta_;
    constantSlipStickCoeffs_.lookup("omega") >> omega_;
    constantSlipStickCoeffs_.lookup("muf") >> muf_;

    return true;
}


// ************************************************************************* //
