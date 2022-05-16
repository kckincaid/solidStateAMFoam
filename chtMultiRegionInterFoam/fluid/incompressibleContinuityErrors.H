{
	volScalarField contErr = fvc::div(phi);

	scalar sumLocalContErr =
	(
		runTime.deltaT().value()*mag(contErr().weightedAverage(mesh.V()).value())
	);

	scalar globalContErr =
	(
		runTime.deltaT().value()*contErr.weightedAverage(mesh.V()).value()
	);

    cumulativeContErr[i] += globalContErr;

    Info<< "time step continuity errors (" << mesh.name() << ")"
        << ": sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr[i]
        << endl;
}
