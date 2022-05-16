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

#include "convectionFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convectionFvPatchScalarField::
convectionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
	kName_("undefined-kName"),
	h_(Zero),
	Tinf_(Zero)
{}


Foam::convectionFvPatchScalarField::
convectionFvPatchScalarField
(
    const convectionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
	kName_(ptf.kName_),
	h_(ptf.h_),
	Tinf_(Zero)
{}


Foam::convectionFvPatchScalarField::
convectionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
	kName_(dict.lookup("k")),
	h_(readScalar(dict.lookup("h"))),
	Tinf_(readScalar(dict.lookup("Tinf")))
{}


Foam::convectionFvPatchScalarField::
convectionFvPatchScalarField
(
    const convectionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
	kName_(ptf.kName_),
	h_(ptf.h_),
	Tinf_(ptf.Tinf_)
{}


Foam::convectionFvPatchScalarField::
convectionFvPatchScalarField
(
    const convectionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
	kName_(ptf.kName_),
	h_(ptf.h_),
	Tinf_(ptf.Tinf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Get thermal conductivity field of fluid
const Foam::fvPatchScalarField&
Foam::convectionFvPatchScalarField::k() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(kName_);

	//Info<< "Flag 2: Looking up mueff name" << endl;
}

void Foam::convectionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Calculate gradient on patch
	tmp<scalarField> f_(1/(1+(k()/(h_*this->patch().deltaCoeffs()))));

	// Calculate field value on patch face
	operator==(this->patchInternalField()*(1-f_) + f_*Tinf_);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::convectionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("kName") << kName_ << token::END_STATEMENT << nl;
	os.writeEntry("h", h_);
	os.writeEntry("Tinf", Tinf_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       convectionFvPatchScalarField
   );
}


// ************************************************************************* //
