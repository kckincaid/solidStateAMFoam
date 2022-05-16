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

#include "constantTemperatureConductionFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantTemperatureConductionFvPatchScalarField::
constantTemperatureConductionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
	L_(Zero),
	Tinf_(Zero)
{}


Foam::constantTemperatureConductionFvPatchScalarField::
constantTemperatureConductionFvPatchScalarField
(
    const constantTemperatureConductionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
	L_(ptf.L_),
	Tinf_(Zero)
{}


Foam::constantTemperatureConductionFvPatchScalarField::
constantTemperatureConductionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
	L_(readScalar(dict.lookup("L"))),
	Tinf_(readScalar(dict.lookup("Tinf")))
{}


Foam::constantTemperatureConductionFvPatchScalarField::
constantTemperatureConductionFvPatchScalarField
(
    const constantTemperatureConductionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
	L_(ptf.L_),
	Tinf_(ptf.Tinf_)
{}


Foam::constantTemperatureConductionFvPatchScalarField::
constantTemperatureConductionFvPatchScalarField
(
    const constantTemperatureConductionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
	L_(ptf.L_),
	Tinf_(ptf.Tinf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constantTemperatureConductionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Calculate fraction expression
	//tmp<scalarField> f_(1/(1+(k()/(h_*this->patch().deltaCoeffs()))));
	tmp<scalarField> f_(1/(1+L_/this->patch().deltaCoeffs()));

	// Calculate field value on patch face
	operator==(this->patchInternalField()*(1-f_) + f_*Tinf_);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::constantTemperatureConductionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
	os.writeEntry("L", L_);
	os.writeEntry("Tinf", Tinf_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       constantTemperatureConductionFvPatchScalarField
   );
}


// ************************************************************************* //
