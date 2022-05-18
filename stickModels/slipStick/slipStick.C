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
    return mu0_*alpha_;
}

// Calculate frictional heating at cell centers
Foam::tmp<Foam::volScalarField>
Foam::stickModels::slipStick::calcQfric() const
{
    return mu0_*alpha_;
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
    mu0_("mu0", dimless, slipStickCoeffs_),
    lambda_("lambda", dimless, slipStickCoeffs_),
    omega_("omega", dimless/dimTime, slipStickCoeffs_),
    etaf_("etaf", dimless, slipStickCoeffs_),
    TC0_("TC0", dimPressure, slipStickCoeffs_),
    TC1_("TC1", dimPressure/dimTemperature, slipStickCoeffs_),
    TC2_("TC2", dimPressure/pow(dimTemperature,2), slipStickCoeffs_),
    TC3_("TC3", dimPressure/pow(dimTemperature,3), slipStickCoeffs_),
    TC4_("TC4", dimPressure/pow(dimTemperature,4), slipStickCoeffs_),

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

    slipStickCoeffs_.lookup("mu0") >> mu0_;
    slipStickCoeffs_.lookup("lambda") >> lambda_;
    slipStickCoeffs_.lookup("omega") >> omega_;
    slipStickCoeffs_.lookup("etaf") >> etaf_;
	slipStickCoeffs_.lookup("TC0") >> TC0_;
	slipStickCoeffs_.lookup("TC1") >> TC1_;
	slipStickCoeffs_.lookup("TC2") >> TC2_;
	slipStickCoeffs_.lookup("TC3") >> TC3_;
	slipStickCoeffs_.lookup("TC4") >> TC4_;

    return true;
}


// ************************************************************************* //
