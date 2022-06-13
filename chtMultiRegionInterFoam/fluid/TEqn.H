{
	// Calculate equivalent strain rate
	epsilonDotEq = Foam::sqrt(2.0/3.0)*mag(symm(fvc::grad(U)));

	// Update effective viscosity
	muEff = thermo.mu();

	// Update heating terms
	friction.correct();

	// Update mixture specific heat (not used in internal calcs)
	cpf = thermo.cp();
	
	// Update mixture conductivity (at cell centers)
	kf = thermo.k();

	// Construct and solve temperature equation
	tmp<fvScalarMatrix> tTEqn
	(
		fvm::ddt(rhoCp, temp)
	  + fvm::div(rhoPhiCp, temp)
	  - fvm::laplacian(kf, temp)
	  //- fvm::Sp(contErr, T) // from compressible interFoam
	 ==
		friction.qvisc()
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
