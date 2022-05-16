// Calculate stress and approximate strain
{	
	//Info << "Beginning strain calculation" << endl;

	// Calculate Zener-Hollomon parameter and flow stress
	dimensionedScalar Q("Q", dimEnergy/dimMoles, thermo.nuModel1().viscosityProperties().lookup("Q"));
	dimensionedScalar R("R", dimEnergy/dimMoles/dimTemperature, thermo.nuModel1().viscosityProperties().lookup("R"));
	dimensionedScalar sigmaR("sigmaR", dimPressure, thermo.nuModel1().viscosityProperties().lookup("sigmaR"));
	dimensionedScalar A("A", dimless/dimTime, thermo.nuModel1().viscosityProperties().lookup("A"));
	dimensionedScalar n("n", dimless, thermo.nuModel1().viscosityProperties().lookup("n"));

	Z = epsilonDotEq*exp(Q/(R*temp));
	sigmaf = sigmaR*asinh(pow(Z/A, 1/n));

	// STRAIN CALCULATION
	// First, calculate incremental strain from current timestep
	volTensorField gradU = fvc::grad(U);
	volTensorField dEpsilon = 0.5*(gradU + gradU.T())*runTime.deltaT();

	// Calculate equivalent plastic strain increment as function of the 2nd invariant of total strain increment
	// Deviatoric component of strain tensor:
	volTensorField dDevEpsilon = dEpsilon - 1.0/3.0*tr(dEpsilon)*I;

	//Info << "Max. dDevEpsilon magnitude: " << max(mag(dDevEpsilon)) << endl;
	// Second invariant of deviator component:
	volScalarField J2s = 0.5*(pow(tr(dDevEpsilon),2) - (dDevEpsilon && dDevEpsilon));

	// Keep only positive values of J2, negatives are rounding error and blow up the sqrt function
	volScalarField J2 = sqrt(pow(J2s,2));

	//Info << "Max. J2 value: " << max(J2).value() << endl;
	//Info << "Min. J2 value: " << min(J2).value() << endl;

	// Equivalent plastic strain increment
	volScalarField dEpsilonEq = sqrt(4.0/3.0*J2);

	Info << "Max. dEpsilonEq value: " << max(dEpsilonEq).value() << endl;

	// Add eq. plastic strain increment to total eq. plastic strain
	epsilonEq += dEpsilonEq;

	Info << "Max. epsilonEq value prior to solving transport eqn: " << max(epsilonEq).value() << endl;

	// Solve transport equation for accumulated strain
	tmp<fvScalarMatrix> tepsilonEqn
	(
		fvm::ddt(rho, epsilonEq)
	  + fvm::div(rhoPhi, epsilonEq)
	  - fvm::laplacian(rho*gammaEpsilon, epsilonEq)
	);

	fvScalarMatrix& epsilonEqn = tepsilonEqn.ref();

	//epsilonEqn.relax();

	epsilonEqn.solve(mesh.solver(epsilonEq.select(finalIter)));

	// Make sure strain field is consistent with BCs
	//epsilonEq.correctBoundaryConditions();

	Info << "Max. epsilonEq value after solving transport eqn: " << max(epsilonEq).value() << endl;

	// Throw away strain values in areas of low alpha for post-processing purposes
	// (don't care about total strain of air, it would be large and screw up rendering)
	//for (int i = 1; i <= epsilonEq.size(); i++)
	//{
	//	if (alpha1[i] < 0.9)
	//	{
	//		epsilonEq[i] = 0;
	//	}
	//}

	Info << "Max. accumulated equivalent strain: " << max(epsilonEq).value() << endl;

	// STRESS CALCULATION
	// First calculate the strain rate
	volTensorField epsilonDot = 0.5*(gradU + gradU.T());
	
	// Keep only deviatoric component of strain rate
	volTensorField epsilonDotDev = dev(epsilonDot);

	// Calculate shear stress tensor using strain rate and viscosity
	volTensorField shearStress = 2.0*muEff*epsilonDotDev;

	// Calculate total stress tensor from shear stress and pressure
	volTensorField sigmaTot = p_rgh*I + shearStress;

	// Calculate deviatoric component of total stress
	volTensorField sigmadev = sigmaTot - 1.0/3.0*tr(sigmaTot)*I;

	// Calculate von Mises stress from deviatoric stress
	sigmavm = Foam::sqrt(2.3/3.0)*Foam::sqrt(sigmadev && sigmadev);

	Info << "Max. von Mises stress: " << max(sigmavm).value() << " Pa" << endl;
}
