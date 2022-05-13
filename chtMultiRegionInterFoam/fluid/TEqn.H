{
	// Calculate frictional heating for boundary conditions
	#include "calculateFrictionalHeating.H"

	// Conversion factor since conductivity fit function is in Celsius
	dimensionedScalar celsConv("celsConv", dimTemperature, -273.0);

	// Calculate effective viscosity
	muEff = thermo.mu();

	// Calculate equivalent strain rate
	epsilonDotEq = Foam::sqrt(2.0/3.0)*mag(symm(fvc::grad(U)));

	// Calculate viscous heating / plastic dissipation term
	qvisc = 3.0*phiVisc*muEff*pow(epsilonDotEq, 2);
	
	// Find fluid conducticity at cell centers using volume average
	kf = thermo.k1()*alpha1 + thermo.k2()*alpha2;

	// Initialize face conductivity outside of if statement
	surfaceScalarField kfh("kf", fvc::interpolate(kf));

	// Recalculate using harmonic interpolation if specified in case properties
	if (interpolationType == "harmonic")
	{
		kfh = 1.0/fvc::interpolate(1/kf);
	}

	// Construct and solve temperature equation
	tmp<fvScalarMatrix> tTEqn
	(
		fvm::ddt(rhoCp, temp)
	  + fvm::div(rhoPhiCp, temp)
	  - fvm::laplacian(kfh, temp)
	  //- fvm::Sp(contErr, T) // from compressible interFoam
	 ==
	  	qvisc
	  + rad.ST(temp)
	  + fvOptions(rhoCp, temp)
	);

	fvScalarMatrix& TEqn = tTEqn.ref();

	TEqn.relax();

	fvOptions.constrain(TEqn);

	TEqn.solve();
    
	fvOptions.correct(temp);

    //thermo.correct();
    rad.correct();

    Info<< "Min/max T:" << min(temp).value() << ' '
        << max(temp).value() << endl;
}
