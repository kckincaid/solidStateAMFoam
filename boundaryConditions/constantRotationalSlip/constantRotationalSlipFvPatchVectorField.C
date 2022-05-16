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

#include "constantRotationalSlipFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantRotationalSlipFvPatchVectorField::
constantRotationalSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    omega_(Zero),
	delta_(Zero),
	Vt_(Zero)
{}


Foam::constantRotationalSlipFvPatchVectorField::
constantRotationalSlipFvPatchVectorField
(
    const constantRotationalSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    omega_(ptf.omega_),
	delta_(ptf.delta_),
	Vt_(ptf.Vt_)
{}


Foam::constantRotationalSlipFvPatchVectorField::
constantRotationalSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    omega_(dict.lookup("omega")),
	delta_(readScalar(dict.lookup("delta"))),
	Vt_(dict.lookup("Vt"))
{}


Foam::constantRotationalSlipFvPatchVectorField::
constantRotationalSlipFvPatchVectorField
(
    const constantRotationalSlipFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    omega_(ptf.omega_),
	delta_(ptf.delta_),
	Vt_(ptf.Vt_)
{}


Foam::constantRotationalSlipFvPatchVectorField::
constantRotationalSlipFvPatchVectorField
(
    const constantRotationalSlipFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    omega_(ptf.omega_),
	delta_(ptf.delta_),
	Vt_(ptf.Vt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constantRotationalSlipFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //const scalar t = this->db().time().timeOutputValue();
    //const scalar axialVelocity = axialVelocity_->value(t);
    //const scalar radialVelocity = radialVelocity_->value(t);
    //const scalar rpm = rpm_->value(t);

    //const vector axisHat = axis_/mag(axis_);
	const vector omegaHat = omega_/mag(omega_);

    //const vectorField r(patch().Cf() - origin_);
	const vectorField r(patch().Cf());

    //const vectorField d(r - (axisHat & r)*axisHat);
	const vectorField d(r - (omegaHat & r)*omegaHat);

    tmp<vectorField> tangVel
    (
        omega_ ^ d
    );

    //operator==(tangVel + axisHat*axialVelocity + radialVelocity*d/mag(d));
	operator==(tangVel*(1-delta_) + Vt_*delta_);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::constantRotationalSlipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("omega", omega_);
	os.writeEntry("delta", delta_);
	os.writeEntry("Vt", Vt_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       constantRotationalSlipFvPatchVectorField
   );
}


// ************************************************************************* //
