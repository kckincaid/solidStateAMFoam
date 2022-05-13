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

#include "incompressibleTwoPhaseThermalMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseThermalMixture, 0);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::incompressibleTwoPhaseThermalMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseThermalMixture::incompressibleTwoPhaseThermalMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rho1_("rho", dimDensity, nuModel1_->viscosityProperties()),
    rho2_("rho", dimDensity, nuModel2_->viscosityProperties()),

	// Constant cubic thermal conductivity values
	k1_("k", dimPower/dimLength/dimTemperature, nuModel1_->viscosityProperties()),
	k2_("k", dimPower/dimLength/dimTemperature, nuModel2_->viscosityProperties()),

	// Temperature-dependent cubic thermal conductivity values
	k10_("k0", dimPower/dimLength/dimTemperature, nuModel1_->viscosityProperties()),
	k11_("k1", dimPower/dimLength/dimTemperature/dimTemperature, nuModel1_->viscosityProperties()),
	k12_("k2", dimPower/dimLength/dimTemperature/dimTemperature/dimTemperature, nuModel1_->viscosityProperties()),
	k13_("k3", dimPower/dimLength/dimTemperature/dimTemperature/dimTemperature/dimTemperature, nuModel1_->viscosityProperties()),

	k20_("k0", dimPower/dimLength/dimTemperature, nuModel2_->viscosityProperties()),
	k21_("k1", dimPower/dimLength/dimTemperature/dimTemperature, nuModel2_->viscosityProperties()),
	k22_("k2", dimPower/dimLength/dimTemperature/dimTemperature/dimTemperature, nuModel2_->viscosityProperties()),
	k23_("k3", dimPower/dimLength/dimTemperature/dimTemperature/dimTemperature/dimTemperature, nuModel2_->viscosityProperties()),

	// Constant specific heats
	cp1_("cp", dimEnergy/dimMass/dimTemperature, nuModel1_->viscosityProperties()),
	cp2_("cp", dimEnergy/dimMass/dimTemperature, nuModel2_->viscosityProperties()),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar(dimViscosity, Zero),
        calculatedFvPatchScalarField::typeName
    )
{
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseThermalMixture::mu() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
          //  limitedAlpha1*rho1_*nuModel1_->nu()
          //+ (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
            // Modified to prevent miscalculation of mu in areas of low alpha1 but large nu1
            limitedAlpha1*limitedAlpha1*rho1_*nuModel1_->nu()
          + (scalar(1) - limitedAlpha1)*(scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
        )
    );
}

Foam::tmp<Foam::scalarField>
Foam::incompressibleTwoPhaseThermalMixture::mu(const label patchI) const
{

    return mu()().boundaryField()[patchI];
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseThermalMixture::muf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muf",
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseThermalMixture::nuf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
            (
                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
              + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
            )/(alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)
        )
    );
}

// Calculate first fluid phase thermal conductivity
Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseThermalMixture::k1() const
{
	// If cubic model, calculate based on temperature
	if ( kModel1_ == "cubic" )
	{
		// Import temperature field
		const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");
		
		return tmp<volScalarField>
    	(
		    new volScalarField
		    (
		        "k1",
				k10_*(temp_/temp_) + k11_*temp_ + k12_*pow(temp, 2) + k13_*pow(temp_,3)
		    )
    	);
	}
	// Otherwise, use constant value
	else
	{
		return tmp<volScalarField>
    	(
		    new volScalarField
		    (
		        "k1",
				k1_*(temp_/temp_)
		    )
    	);
	}
}

// Calculate second fluid phase thermal conductivity
Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseThermalMixture::k2() const
{
	// If cubic model, calculate based on temperature
	if ( kModel2_ == "cubic" )
	{
		// Import temperature field
		const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");
		
		return tmp<volScalarField>
    	(
		    new volScalarField
		    (
		        "k2",
				k20_*(temp_/temp_) + k21_*temp_ + k22_*pow(temp, 2) + k23_*pow(temp_,3)
		    )
    	);
	}
	// Otherwise, use constant value
	else
	{
		return tmp<volScalarField>
    	(
		    new volScalarField
		    (
		        "k2",
				k2_*(temp_/temp_)
		    )
    	);
	}
}

// Calculate fluid mixture thermal conductivity
Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseThermalMixture::k() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "k",
            limitedAlpha1*k1() + (scalar(1) - limitedAlpha1)*k2()
        )
    );
}

// Calculate fluid mixture specific heat
Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseThermalMixture::cp() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cp",
            limitedAlpha1*cp1_ + (scalar(1) - limitedAlpha1)*cp2_
        )
    );
}

// Calculate fluid mixture density
Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseThermalMixture::rho() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rho",
            limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_
        )
    );
}

bool Foam::incompressibleTwoPhaseThermalMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            nuModel1_().read
            (
                subDict(phase1Name_ == "1" ? "phase1": phase1Name_)
            )
         && nuModel2_().read
            (
                subDict(phase2Name_ == "2" ? "phase2": phase2Name_)
            )
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;

			// Read thermal conductivity types
			nuModel1_->viscosityProperties().lookup("kModel") >> kModel1_;
			nuModel2_->viscosityProperties().lookup("kModel") >> kModel2_;

			// Look up appropriate constants based on model type (constant or cubic)
			if ( kModel1_ == "cubic" )
			{
				nuModel1_->viscosityProperties().lookup("k0") >> k10_;
				nuModel1_->viscosityProperties().lookup("k1") >> k11_;
				nuModel1_->viscosityProperties().lookup("k2") >> k12_;
				nuModel1_->viscosityProperties().lookup("k3") >> k13_;
			}
			else
			{
				nuModel1_->viscosityProperties().lookup("k") >> k1_;
			}

			// Repeat for second fluid phase
			if ( kModel2_ == "cubic" )
			{
				nuModel2_->viscosityProperties().lookup("k0") >> k20_;
				nuModel2_->viscosityProperties().lookup("k1") >> k21_;
				nuModel2_->viscosityProperties().lookup("k2") >> k22_;
				nuModel2_->viscosityProperties().lookup("k3") >> k23_;
			}
			else
			{
				nuModel2_->viscosityProperties().lookup("k") >> k2_;
			}

			// Read in specific heats for each phase
			nuModel1_->viscosityProperties().lookup("cp") >> cp1_;
			nuModel2_->viscosityProperties().lookup("cp") >> cp2_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
