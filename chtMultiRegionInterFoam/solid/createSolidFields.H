    // Initialise solid field pointer lists
    PtrList<radiation::radiationModel> radiations(solidRegions.size());
    PtrList<fv::options> solidHeatSources(solidRegions.size());
    PtrList<volScalarField> TSolid(solidRegions.size());
	PtrList<IOdictionary> solidProps(solidRegions.size());
	PtrList<dimensionedScalar> rhoSolid(solidRegions.size());
	PtrList<dimensionedScalar> cpSolid(solidRegions.size());
	PtrList<volScalarField> kSolid(solidRegions.size());

	// Solid velocity value and field
	PtrList<volVectorField> USolid(solidRegions.size());
	PtrList<surfaceScalarField> phiSolid(solidRegions.size());

	// Temperature-dependent thermal conductivity values
	PtrList<dimensionedScalar> k0Solid(solidRegions.size());
	PtrList<dimensionedScalar> k1Solid(solidRegions.size());
	PtrList<dimensionedScalar> k2Solid(solidRegions.size());
	PtrList<dimensionedScalar> k3Solid(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding fvOptions\n" << endl;
        solidHeatSources.set
        (
            i,
            new fv::options(solidRegions[i])
        );

		Info<< "    Adding to TSolid\n" << endl;
        TSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "temp",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

		Info<< "    Adding to radiations\n" << endl;
        radiations.set(i, radiation::radiationModel::New(TSolid[i]));

		// Read solidProperties file
		Info << "	reading solidProperties\n" << endl;

		solidProps.set
		(
			i,
			new IOdictionary
			(
				IOobject
				(
					"solidProperties",
					runTime.constant(),
					solidRegions[i],
					IOobject::MUST_READ_IF_MODIFIED,
					IOobject::NO_WRITE
				)
			)
		);

		// Solid density
		rhoSolid.set
		(
			i,
			new dimensionedScalar
			(
				"rho",
				dimDensity,
				solidProps[i].lookup("rho")
			)
		);

		// Solid specific heat
		cpSolid.set
		(
			i,
			new dimensionedScalar
			(
				"cp",
				dimEnergy/dimMass/dimTemperature,
				solidProps[i].lookup("cp")
			)
		);

		// Solid thermal conductivity (scalar)
		k0Solid.set
		(
			i,
			new dimensionedScalar
			(
				"k0",
				dimensionSet(1, 1, -3, -1, 0, 0, 0),
				solidProps[i].lookup("k0")
			)
		);

		k1Solid.set
		(
			i,
			new dimensionedScalar
			(
				"k1",
				dimensionSet(1, 1, -3, -2, 0, 0, 0),
				solidProps[i].lookup("k1")
			)
		);

		k2Solid.set
		(
			i,
			new dimensionedScalar
			(
				"k2",
				dimensionSet(1, 1, -3, -3, 0, 0, 0),
				solidProps[i].lookup("k2")
			)
		);

		k3Solid.set
		(
			i,
			new dimensionedScalar
			(
				"k3",
				dimensionSet(1, 1, -3, -4, 0, 0, 0),
				solidProps[i].lookup("k3")
			)
		);

		// Solid thermal conductivity (volScalarField)
        Info<< "    Adding to kSolid\n" << endl;
        kSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "k",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
				k0Solid[i]
            )
        );
		
		USolid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
				Vt //VSolid[i]
            )
        );

		Info<< "    Adding to flux phiSolid\n" << endl;
		phiSolid.set
		(
		    i,
		    new surfaceScalarField
		    (
		        IOobject
		        (
		            "phi",
		            runTime.timeName(),
		            solidRegions[i],
		            IOobject::READ_IF_PRESENT,
		            IOobject::NO_WRITE
		        ),
		        linearInterpolate(USolid[i]) & solidRegions[i].Sf()
		    )
		);
    }
