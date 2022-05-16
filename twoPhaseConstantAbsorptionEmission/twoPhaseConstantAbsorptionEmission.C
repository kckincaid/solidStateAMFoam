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

#include "twoPhaseConstantAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(twoPhaseConstantAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            twoPhaseConstantAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::twoPhaseConstantAbsorptionEmission::twoPhaseConstantAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    a1_(coeffsDict_.lookup("absorptivity1")),
	a2_(coeffsDict_.lookup("absorptivity2")),
    e1_(coeffsDict_.lookup("emissivity1")),
	e2_(coeffsDict_.lookup("emissivity2")),
    E1_(coeffsDict_.lookup("E1")),
	E2_(coeffsDict_.lookup("E2")),
	alphaField_(coeffsDict_.lookup("alpha"))
	//alpha_(mesh.lookupObject<volScalarField>("alpha1"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::twoPhaseConstantAbsorptionEmission::~twoPhaseConstantAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::twoPhaseConstantAbsorptionEmission::aCont(const label bandI) const
{
	// Import phase fraction
	const volScalarField& alpha_ = mesh_.lookupObject<volScalarField>(alphaField_);

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            a1_
        )
    );

	ta = alpha_*a1_ + (scalar(1) - alpha_)*a2_;

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::twoPhaseConstantAbsorptionEmission::eCont(const label bandI) const
{
	// Import phase fraction
	const volScalarField& alpha_ = mesh_.lookupObject<volScalarField>(alphaField_);

    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            e1_
        )
    );

	te = alpha_*e1_ + (scalar(1) - alpha_)*e2_;

    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::twoPhaseConstantAbsorptionEmission::ECont(const label bandI) const
{
	// Import phase fraction
	const volScalarField& alpha_ = mesh_.lookupObject<volScalarField>(alphaField_);

    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            E1_
        )
    );

	tE = alpha_*E1_ + (scalar(1) - alpha_)*E2_;

    return tE;
}


// ************************************************************************* //
