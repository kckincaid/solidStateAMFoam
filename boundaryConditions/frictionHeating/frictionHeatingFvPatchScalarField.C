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

#include "frictionHeatingFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frictionHeatingFvPatchScalarField::
frictionHeatingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
	qfricName_("undefined-qfricName"),
	kName_("undefined-kName"),
	f_(Zero)
{}


Foam::frictionHeatingFvPatchScalarField::
frictionHeatingFvPatchScalarField
(
    const frictionHeatingFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
	qfricName_(ptf.qfricName_),
	kName_(ptf.kName_),
	f_(ptf.f_)
{}


Foam::frictionHeatingFvPatchScalarField::
frictionHeatingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
	qfricName_(dict.lookup("qfric")),
	kName_(dict.lookup("k")),
	f_(readScalar(dict.lookup("f")))
{}


Foam::frictionHeatingFvPatchScalarField::
frictionHeatingFvPatchScalarField
(
    const frictionHeatingFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
	qfricName_(ptf.qfricName_),
	kName_(ptf.kName_),
	f_(ptf.f_)
{}


Foam::frictionHeatingFvPatchScalarField::
frictionHeatingFvPatchScalarField
(
    const frictionHeatingFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
	qfricName_(ptf.qfricName_),
	kName_(ptf.kName_),
	f_(ptf.f_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Get pressure field from corresponding pName entry
const Foam::fvPatchScalarField&
Foam::frictionHeatingFvPatchScalarField::friction() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(qfricName_);

	//Info<< "Flag 1: Looking up p name" << endl;
}

// Get effective viscosity field from corresponding mueffName Entry
const Foam::fvPatchScalarField&
Foam::frictionHeatingFvPatchScalarField::k() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(kName_);

	//Info<< "Flag 2: Looking up mueff name" << endl;
}

void Foam::frictionHeatingFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Calculate gradient on patch
	tmp<scalarField> grad_(f_*friction()/k());

	// Calculate field value on patch face
	operator==(this->patchInternalField() + grad_/this->patch().deltaCoeffs());

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::frictionHeatingFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
	os.writeKeyword("qfricName")<< qfricName_ << token::END_STATEMENT << nl;
    os.writeKeyword("kName") << kName_ << token::END_STATEMENT << nl;
	os.writeEntry("f", f_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       frictionHeatingFvPatchScalarField
   );
}


// ************************************************************************* //
