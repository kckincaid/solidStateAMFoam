/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "incompressibleFrictionModel.H"
#include "stickModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleFrictionModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleFrictionModel::incompressibleFrictionModel
(
    const volVectorField& U,
	const volVectorField& Ut,
    const volScalarField& alpha,
	const volScalarField& p
)
:
    IOdictionary
    (
        IOobject
        (
            "frictionProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    stickModelPtr_(stickModel::New("qvisc", *this, U, Ut, alpha, p))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompressibleFrictionModel::~incompressibleFrictionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleFrictionModel::qvisc() const
{
    return stickModelPtr_->qvisc();
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleFrictionModel::qfric() const
{
    return stickModelPtr_->qfric();
}


void Foam::incompressibleFrictionModel::correct()
{
    stickModelPtr_->correct();
}


bool Foam::incompressibleFrictionModel::read()
{
    if (regIOobject::read())
    {
        return stickModelPtr_->read(*this);
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
