/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "NortonHoff.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(NortonHoff, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        NortonHoff,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::NortonHoff::calcNu() const
{
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
			// Norton-Hoff viscosity calculation
			(mu0_/rho_)*pow
						(
							sqrt(3.0)*max(strainRate(),strainRateEffMin_)/strainRateUnits_,
							(m_-1)
						)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::NortonHoff::NortonHoff
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    NortonHoffCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    mu0_("mu0", dimViscosity*dimDensity, NortonHoffCoeffs_),
	m_("m", dimless, NortonHoffCoeffs_),
	rho_("rho", dimDensity, NortonHoffCoeffs_),
    nuMin_("nuMin", dimViscosity, NortonHoffCoeffs_),
    nuMax_("nuMax", dimViscosity, NortonHoffCoeffs_),
	strainRateEffMin_("strainRateEffMin", dimless/dimTime, NortonHoffCoeffs_),
	strainRateUnits_("strainRateUnits", dimless/dimTime, 1),

    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::NortonHoff::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    NortonHoffCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    NortonHoffCoeffs_.lookup("mu0") >> mu0_;
    NortonHoffCoeffs_.lookup("m") >> m_;
	NortonHoffCoeffs_.lookup("rho") >> rho_;
    NortonHoffCoeffs_.lookup("nuMin") >> nuMin_;
    NortonHoffCoeffs_.lookup("nuMax") >> nuMax_;
	NortonHoffCoeffs_.lookup("strainRateEffMin") >> strainRateEffMin_;

    return true;
}


// ************************************************************************* //
