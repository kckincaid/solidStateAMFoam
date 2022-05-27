if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
	dimensionedScalar celsConv("celsConv", dimTemperature, -273.0);

	// Calculate temp-dependent thermal conductivity
	k = k0s + k1s*(temp+celsConv) + k2s*pow(temp+celsConv,2) + k3s*pow(temp+celsConv,3);
	
	// Limit thermal conductivity to prevent instability
	forAll(temp.internalField(), celli)
	{
		// If temp is above maximum temperature, limit conductivity to value at Tmax
		if ( temp[celli] > Tmax.value() )
		{
			k[celli] = (k0s*(Tmax/Tmax) + k1s*(Tmax+celsConv) + k2s*pow(Tmax+celsConv,2) + k3s*pow(Tmax+celsConv,3)).value();
		}
		
		// If temp is below minumum temperature, limit conductivity to value at Tmin
		if ( temp[celli] < Tmin.value() )
		{
			k[celli] = (k0s*(Tmin/Tmin) + k1s*(Tmin+celsConv) + k2s*pow(Tmin+celsConv,2) + k3s*pow(Tmin+celsConv,3)).value();
		}
	}

	// Initialize face conductivity outside of if statement
	surfaceScalarField ksf("kf", fvc::interpolate(k));

	// Recalculate using harmonic interpolation if specified in case properties
	if (interpolationType == "harmonic")
	{
		ksf = 1.0/fvc::interpolate(1/k);
	}

	// Solve heat equation in solid
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho*cp, temp)
	  + fvm::div(rho*cp*phi, temp)
	  - fvm::laplacian(ksf, temp)
      ==
        rho*cp*fvOptions(temp)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(temp);

    //thermo.correct();

    Info<< "Min/max T:" << min(temp).value() << ' '
        << max(temp).value() << endl;
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
