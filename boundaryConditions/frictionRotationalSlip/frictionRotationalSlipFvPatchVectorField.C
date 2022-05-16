/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "frictionRotationalSlipFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frictionRotationalSlipFvPatchVectorField::
frictionRotationalSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
	pName_("undefined-pName"),
	mueffName_("undefined-mueffName"),
	muf_(Zero),
	omega_(Zero),
	Vt_(Zero)
{}


Foam::frictionRotationalSlipFvPatchVectorField::
frictionRotationalSlipFvPatchVectorField
(
    const frictionRotationalSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
	pName_(ptf.pName_),
	mueffName_(ptf.mueffName_),
	muf_(ptf.muf_),
    omega_(ptf.omega_),
	Vt_(ptf.Vt_)
{}


Foam::frictionRotationalSlipFvPatchVectorField::
frictionRotationalSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
	pName_(dict.lookup("p")),
	mueffName_(dict.lookup("mueff")),
	muf_(readScalar(dict.lookup("muf"))),
    omega_(dict.lookup("omega")),
	Vt_(dict.lookup("Vt"))
{}


Foam::frictionRotationalSlipFvPatchVectorField::
frictionRotationalSlipFvPatchVectorField
(
    const frictionRotationalSlipFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
	pName_(ptf.pName_),
	mueffName_(ptf.mueffName_),
	muf_(ptf.muf_),
    omega_(ptf.omega_),
	Vt_(ptf.Vt_)
{}


Foam::frictionRotationalSlipFvPatchVectorField::
frictionRotationalSlipFvPatchVectorField
(
    const frictionRotationalSlipFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
	pName_(ptf.pName_),
	mueffName_(ptf.mueffName_),
	muf_(ptf.muf_),
    omega_(ptf.omega_),
	Vt_(ptf.Vt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Get pressure field from corresponding pName entry
const Foam::fvPatchScalarField&
Foam::frictionRotationalSlipFvPatchVectorField::pressure() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(pName_);

	//Info<< "Flag 1: Looking up p name" << endl;
}

// Get effective viscosity field from corresponding mueffName Entry
const Foam::fvPatchScalarField&
Foam::frictionRotationalSlipFvPatchVectorField::mueff() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(mueffName_);

	//Info<< "Flag 2: Looking up mueff name" << endl;
}

void Foam::frictionRotationalSlipFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //const vector axisHat = axis_/mag(axis_);
	const vector omegaHat = omega_/mag(omega_);

    //const vectorField r(patch().Cf() - origin_);
	const vectorField r(patch().Cf());

    //const vectorField d(r - (axisHat & r)*axisHat);
	const vectorField d(r - (omegaHat & r)*omegaHat);

	//const scalarField R(mag(d));

	//const scalarField delta(1.0 - exp(-mag(omega_)*R/(delta0_*omega0_*R0_)));

	// Calculate tangential velocity field on patch
    	//tmp<vectorField> tangVel(omega_ ^ d);
	const vectorField tangVel(omega_ ^ d);

	// Calculate tangential velocity unit vector on patch
	//tmp<vectorField> tangVelHat(tangVel/mag(tangVel));
	const vectorField tangVelHat(tangVel/mag(tangVel));

	// Calculate velocity gradient
	//tmp<vectorField> grad_(tangVelHat*muf_*pressure()/mueff());
	const vectorField grad_(tangVelHat*muf_*pressure()/max(mueff(),VSMALL));

    //operator==(tangVel + axisHat*axialVelocity + radialVelocity*d/mag(d));
	//operator==(tangVel*(1-delta) + Vt_*delta);

	// Calculate field value on patch face
	operator==(this->patchInternalField() + grad_/this->patch().deltaCoeffs());

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::frictionRotationalSlipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
	os.writeKeyword("pName")<< pName_ << token::END_STATEMENT << nl;
    os.writeKeyword("mueffName") << mueffName_ << token::END_STATEMENT << nl;
	os.writeEntry("muf", muf_);
	os.writeEntry("omega", omega_);
	os.writeEntry("Vt", Vt_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       frictionRotationalSlipFvPatchVectorField
   );
}


// ************************************************************************* //
