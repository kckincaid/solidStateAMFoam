{
	// Calculate frictional heating for boundary conditions
	//#include "calculateFrictionalHeating.H"

	// Calculate equivalent strain rate
	epsilonDotEq = Foam::sqrt(2.0/3.0)*mag(symm(fvc::grad(U)));

	// Update effective viscosity
	muEff = thermo.mu();

	// Update heating terms
	friction.correct();
	qfric = friction.qfric();
	qvisc = friction.qvisc();

	// Update qfric on boundary faces
	qfric.correctBoundaryConditions();

	// Update mixture specific heat (not used in internal calcs)
	cpf = thermo.cp();
	
	// Update mixture conductivity (at cell centers)
	kf = thermo.k();

	// Initialize face conductivity with arithmetic interpolation
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
		//friction.qvisc()
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
