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

	// Constant densities
    rho1_("rho", dimDensity, nuModel1_->viscosityProperties()),
    rho2_("rho", dimDensity, nuModel2_->viscosityProperties()),

	// Constant specific heats
	cp1_("cp", dimEnergy/dimMass/dimTemperature, nuModel1_->viscosityProperties()),
	cp2_("cp", dimEnergy/dimMass/dimTemperature, nuModel2_->viscosityProperties()),

	// Thermal conductivity models and units
	kModel1_(nuModel1_->viscosityProperties().lookupOrDefault<word>("kModel","constant")),
	kModel2_(nuModel2_->viscosityProperties().lookupOrDefault<word>("kModel","constant")),

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
	// Import temperature field
	const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");

	// If cubic model, calculate based on temperature
	if ( kModel1_ == "cubic" )
	{
		// Look up necessary constants
		dimensionedScalar k10_("k0", dimPower/dimLength/dimTemperature, nuModel1_->viscosityProperties().lookup("k0"));
		dimensionedScalar k11_("k1", dimPower/dimLength/dimTemperature/dimTemperature, nuModel1_->viscosityProperties().lookup("k1"));
		dimensionedScalar k12_("k2", dimPower/dimLength/dimTemperature/dimTemperature/dimTemperature, nuModel1_->viscosityProperties().lookup("k2"));
		dimensionedScalar k13_("k3", dimPower/dimLength/dimTemperature/dimTemperature/dimTemperature/dimTemperature, nuModel1_->viscosityProperties().lookup("k3"));
		dimensionedScalar Tmax1_("Tmax", dimTemperature, nuModel1_->viscosityProperties().lookup("Tmax"));	
		dimensionedScalar Tmin1_("Tmin", dimTemperature, nuModel1_->viscosityProperties().lookup("Tmin"));

		// Look up units of fit function (Celsius or Kelvin)
		word units1_(nuModel1_->viscosityProperties().lookupOrDefault<word>("units","Kelvin"));

		// Declare conversion variable prior to setting value
		dimensionedScalar celsConv_("celsConv", dimTemperature, -273.0);

		// Set conversion to zero if cubic fit is already in Kelvin
		if ( units1_ == "Kelvin" )
		{
			celsConv_ = 0.0;
		}
		else if ( units1_ != "Celsius" )
		{
            FatalErrorInFunction
                << "Unknown conductivity function units " << units1_
				<< " for phase 1. Choose Kelvin or Celsius."
                << exit(FatalError);
		}
	
		// Convert temp field if needed and limit to prevent inaccuracy at high/low temps
		const volScalarField limitedTemp_
		(
		    min(max(temp_ + celsConv_, Tmin1_), Tmax1_)
		);
		
		// Calculate and return temp-dependent thermal conductivity
		return tmp<volScalarField>
    	(
		    new volScalarField
		    (
		        "k1",
				k10_*(limitedTemp_/limitedTemp_) + k11_*limitedTemp_ + k12_*pow(limitedTemp_, 2) + k13_*pow(limitedTemp_,3)
		    )
    	);
	}
	// If not cubic model, use constant value
	else
	{
		// Look up constant value from dictionary
		dimensionedScalar k1_("k", dimPower/dimLength/dimTemperature, nuModel1_->viscosityProperties().lookup("k"));

		// Return a volScalarField with specified constant value
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
	// Import temperature field
	const volScalarField& temp_ = U_.mesh().lookupObject<volScalarField>("temp");

	// If cubic model, calculate based on temperature
	if ( kModel2_ == "cubic" )
	{
		dimensionedScalar k20_("k0", dimPower/dimLength/dimTemperature, nuModel2_->viscosityProperties().lookup("k0"));
		dimensionedScalar k21_("k1", dimPower/dimLength/dimTemperature/dimTemperature, nuModel2_->viscosityProperties().lookup("k1"));
		dimensionedScalar k22_("k2", dimPower/dimLength/dimTemperature/dimTemperature/dimTemperature, nuModel2_->viscosityProperties().lookup("k2"));
		dimensionedScalar k23_("k3", dimPower/dimLength/dimTemperature/dimTemperature/dimTemperature/dimTemperature, nuModel2_->viscosityProperties().lookup("k3"));
		dimensionedScalar Tmax2_("Tmax", dimTemperature, nuModel2_->viscosityProperties().lookup("Tmax"));	
		dimensionedScalar Tmin2_("Tmin", dimTemperature, nuModel2_->viscosityProperties().lookup("Tmin"));

		// Look up units of fit function (Celsius or Kelvin)
		word units2_(nuModel2_->viscosityProperties().lookupOrDefault<word>("units","Kelvin"));

		// Declare conversion variable prior to setting value
		dimensionedScalar celsConv_("celsConv", dimTemperature, -273.0);

		// Set conversion to zero if cubic fit is already in Kelvin
		if ( units2_ == "Kelvin" )
		{
			celsConv_ = 0.0;
		}
		else if ( units2_ != "Celsius" )
		{
            FatalErrorInFunction
                << "Unknown conductivity function units " << units2_
				<< " for phase 2. Choose Kelvin or Celsius."
                << exit(FatalError);
		}

		// Create limited temperature field to prevent instabilities and convert to Celsius
		const volScalarField limitedTemp_
		(
		    min(max(temp_ + celsConv_, Tmin2_), Tmax2_)
		);
		
		// Calculate and return temp-dependent thermal conductivity
		return tmp<volScalarField>
    	(
		    new volScalarField
		    (
		        "k2",
				k20_*(limitedTemp_/limitedTemp_) + k21_*limitedTemp_ + k22_*pow(limitedTemp_, 2) + k23_*pow(limitedTemp_,3)
		    )
    	);
	}
	// If not a cubic model, use constant value
	else
	{
		// Look up constant value from dictionary
		dimensionedScalar k2_("k", dimPower/dimLength/dimTemperature, nuModel2_->viscosityProperties().lookup("k"));

		// Return a volScalarField with specified constant value
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
