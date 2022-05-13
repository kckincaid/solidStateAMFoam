PtrList<IOobject> alphaPhi10HeaderList(fluidRegions.size());
PtrList<bool> alphaRestartList(fluidRegions.size());
PtrList<surfaceScalarField> alphaPhi10List(fluidRegions.size());
//PtrList<surfaceScalarField> alphaPhi1Corr0List(fluidRegions.size());

Info << "	Creating alpha fluxes\n" << endl;

forAll(fluidRegions, i)
{
	alphaPhi10HeaderList.set
	(
		i,
		new IOobject
		(
			IOobject
			(
				IOobject::groupName("alphaPhi0", alpha1l[i].group()),
				runTime.timeName(),
				fluidRegions[i],
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		)
	);

	//IOobject alphaPhi10Header
	//(
	//    IOobject::groupName("alphaPhi0", alpha1.group()),
	//    runTime.timeName(),
	//    mesh,
	//    IOobject::READ_IF_PRESENT,
	//    IOobject::AUTO_WRITE
	//);

	

	IOobject alphaRestartTmp = alphaPhi10HeaderList[i];

	const bool alphaRestartListTmp = alphaRestartTmp.typeHeaderOk<surfaceScalarField>(true);

	alphaRestartList.set
	(
		i,
		new bool(alphaRestartListTmp)
	);
	//const bool alphaRestartList[i] = alphaRestartListTmp;
		//alphaRestartTmp.typeHeaderOk<surfaceScalarField>(true);

	//const bool alphaRestart =
		//alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true);

	Info << "	Alpha flux flag 1\n" << endl;

	if (alphaRestartList[i])
	{
		Info << "Restarting alpha" << endl;
	}

	// MULES flux from previous time-step

	alphaPhi10List.set
	(
		i,
		new surfaceScalarField
		(
			"alphaPhi10Header",
			phiFluid[i]*fvc::interpolate(alpha1l[i])
		)
	);

	Info << "	Alpha flux flag 2\n" << endl;

	//surfaceScalarField alphaPhi10
	//(
	//	alphaPhi10Header,
	//	phi*fvc::interpolate(alpha1)
	//);

	tmp<surfaceScalarField> talphaPhi1Corr0;

	// MULES Correction
	//alphaPhi1Corr0List.set
	//(
		//i,
		//new surfaceScalarField //tmp<surfaceScalarField> doesn't work
		//(
			//talphaPhi1Corr0()
			//"talphaPhi1Corr0",
			//0 // needed this value to initialize, unsure if correct
		//)
	//);

	Info << "	Done creating alpha fluxes\n" << endl;

	//tmp<surfaceScalarField> talphaPhi1Corr0;
}
