// Initialise fluid field pointer lists

// Transport model
PtrList<immiscibleIncompressibleTwoPhaseThermalMixture> twoPhaseFluid(fluidRegions.size());

// Flow variables
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());

// Gravity terms for pressure calculations
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());

// Weld velocity field for modified reference frame
PtrList<dimensionedVector> VtFluid(fluidRegions.size());
PtrList<volVectorField> UtFluid(fluidRegions.size());

// Calculate effective viscosity, strain rate, and viscous heating
PtrList<IOdictionary> transportPropertiesl(fluidRegions.size());
PtrList<dimensionedScalar> phiViscl(fluidRegions.size());
PtrList<volScalarField> muEffl(fluidRegions.size());
PtrList<volScalarField> epsilonDotEql(fluidRegions.size());
PtrList<volScalarField> qviscl(fluidRegions.size());

// Fluid average physical properties
PtrList<volScalarField> alpha1l(fluidRegions.size());
PtrList<volScalarField> alpha2l(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size()); // mass flux using rho & phi
PtrList<volScalarField> rhoCpFluid(fluidRegions.size()); // need to store rho*cp for ddt(rho*cp,T)
PtrList<surfaceScalarField> rhoPhiCpFluid(fluidRegions.size()); // energy flux
PtrList<volScalarField> kFluid(fluidRegions.size());
PtrList<volScalarField> cpFluid(fluidRegions.size());

// Phase individual properties
PtrList<dimensionedScalar> rhoFluid1(fluidRegions.size());
PtrList<dimensionedScalar> rhoFluid2(fluidRegions.size());
PtrList<dimensionedScalar> cpFluid1(fluidRegions.size());
PtrList<dimensionedScalar> cpFluid2(fluidRegions.size());
PtrList<volScalarField> kFluid1(fluidRegions.size());
PtrList<volScalarField> kFluid2(fluidRegions.size());

// Cubic thermal conductivity coefficients
PtrList<dimensionedScalar> k0Fluid1(fluidRegions.size());
PtrList<dimensionedScalar> k1Fluid1(fluidRegions.size());
PtrList<dimensionedScalar> k2Fluid1(fluidRegions.size());
PtrList<dimensionedScalar> k3Fluid1(fluidRegions.size());
PtrList<dimensionedScalar> k0Fluid2(fluidRegions.size());
PtrList<dimensionedScalar> k1Fluid2(fluidRegions.size());
PtrList<dimensionedScalar> k2Fluid2(fluidRegions.size());
PtrList<dimensionedScalar> k3Fluid2(fluidRegions.size());

// Quartic temp-dependent max shear stress coeffs
// Correlation is for temp (K) to max shear stress (MPa)
PtrList<dimensionedScalar> TC0Fluid(fluidRegions.size());
PtrList<dimensionedScalar> TC1Fluid(fluidRegions.size());
PtrList<dimensionedScalar> TC2Fluid(fluidRegions.size());
PtrList<dimensionedScalar> TC3Fluid(fluidRegions.size());
PtrList<dimensionedScalar> TC4Fluid(fluidRegions.size());
PtrList<dimensionedScalar> TminFluid(fluidRegions.size());
PtrList<dimensionedScalar> TmaxFluid(fluidRegions.size());
PtrList<dimensionedScalar> tauminFluid(fluidRegions.size());
PtrList<dimensionedScalar> taumaxFluid(fluidRegions.size());

// Stress and strain-related variables
List<bool> trackStressStrainl(fluidRegions.size());
PtrList<volScalarField> ZFluid(fluidRegions.size());
PtrList<volScalarField> sigmafFluid(fluidRegions.size());
PtrList<volScalarField> epsilonEqFluid(fluidRegions.size());
PtrList<dimensionedScalar> gammaEpsilonFluid(fluidRegions.size());
PtrList<volScalarField> sigmavmFluid(fluidRegions.size());

// Particle tracking variables
List<bool> particleTrackingl(fluidRegions.size());
PtrList<volScalarField> betaFluid(fluidRegions.size());
PtrList<dimensionedScalar> gammaBetaFluid(fluidRegions.size());

// Celsius conversion
//PtrList<dimensionedScalar> celsConvl(fluidRegions.size());

// Other necessary lists
PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
PtrList<radiation::radiationModel> radiationFluid(fluidRegions.size());
List<label> pRefCellFluid(fluidRegions.size(),0);
List<scalar> pRefValueFluid(fluidRegions.size(),0.0);

// Frictional heat generation terms
List<word> slipTypel(fluidRegions.size());					// Type of condition to apply at slip faces
PtrList<dimensionedScalar> etafl(fluidRegions.size());		// Mechanical efficiency of frictional heating
PtrList<dimensionedScalar> omegal(fluidRegions.size());		// Angular velocity
PtrList<dimensionedScalar> omega0l(fluidRegions.size());	// Model constant
PtrList<dimensionedScalar> delta0l(fluidRegions.size());	// Model constant
PtrList<dimensionedScalar> Rsl(fluidRegions.size());		// Tool shoulder diameter
PtrList<dimensionedScalar> mu0l(fluidRegions.size());		// Model constant
PtrList<dimensionedScalar> lambdal(fluidRegions.size());	// Model constant
PtrList<volScalarField> qfricl(fluidRegions.size());		// Frictional heat field used in temperature BC

// Dynamic Mesh control variables
List<bool> correctPhiList(fluidRegions.size());
List<bool> checkMeshCourantNoList(fluidRegions.size());
List<bool> moveMeshOuterCorrectorsList(fluidRegions.size());

List<bool> frozenFlowFluid(fluidRegions.size(), false);

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to velocity UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to flux phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            linearInterpolate(UFluid[i]) & fluidRegions[i].Sf()
        )
    );
    Info<< "    Adding to temperature field TFluid\n" << endl;
    TFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "temp",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

	Info<< "    Adding to effective viscosity muEff\n" << endl;
    muEffl.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "muEff",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
			fluidRegions[i],
			dimensionedScalar("muEff", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
        )
    );

	Info<< "    Adding to equivalent strain rate epsilonDotEq\n" << endl;
    epsilonDotEql.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "epsilonDotEq",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
			fluidRegions[i],
			dimensionedScalar("epsilonDotEq", dimless/dimTime, 0.0)
        )
    );

    Info<< "    Adding to viscous heating qvisc\n" << endl;
    qviscl.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "qvisc",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
			fluidRegions[i],
			dimensionedScalar("qvisc", dimDensity*dimViscosity/dimTime/dimTime, 0.0)
        )
    );

	// Create two phase transport model for each fluid region
	Info<< "	Creating two phase transport model\n" << endl;
	twoPhaseFluid.set
	(
		i,
		new immiscibleIncompressibleTwoPhaseThermalMixture(UFluid[i], phiFluid[i])
	);

	Info<< "	Adding to phase fraction alpha\n" << endl;

	// Populate volume fraction fields
	alpha1l.set
	(
		i,
		new volScalarField
		(
			"alpha1",
			twoPhaseFluid[i].alpha1()
		)
	);

	alpha2l.set
	(
		i,
		new volScalarField
		(
			"alpha2",
			twoPhaseFluid[i].alpha2()
		)
	);

	rhoFluid1.set
	(
		i,
		new dimensionedScalar
		(
			"rho1",
			twoPhaseFluid[i].rho1()
		)
	);

	rhoFluid2.set
	(
		i,
		new dimensionedScalar
		(
			"rho2",
			twoPhaseFluid[i].rho2()
		)
	);

	cpFluid1.set
	(
		i,
		new dimensionedScalar
		(
			"cp1",
			twoPhaseFluid[i].cp1()
		)
	);

	cpFluid2.set
	(
		i,
		new dimensionedScalar
		(
			"cp2",
			twoPhaseFluid[i].cp2()
		)
	);

	k0Fluid1.set
	(
		i,
		new dimensionedScalar
		(
			"k0_1",
			dimensionSet(1, 1, -3, -1, 0, 0, 0),
			twoPhaseFluid[i].nuModel1().viscosityProperties().lookup("k0")
		)
	);

	k1Fluid1.set
	(
		i,
		new dimensionedScalar
		(
			"k1_1",
			dimensionSet(1, 1, -3, -2, 0, 0, 0),
			twoPhaseFluid[i].nuModel1().viscosityProperties().lookup("k1")
		)
	);

	k2Fluid1.set
	(
		i,
		new dimensionedScalar
		(
			"k2_1",
			dimensionSet(1, 1, -3, -3, 0, 0, 0),
			twoPhaseFluid[i].nuModel1().viscosityProperties().lookup("k2")
		)
	);

	k3Fluid1.set
	(
		i,
		new dimensionedScalar
		(
			"k3_1",
			dimensionSet(1, 1, -3, -4, 0, 0, 0),
			twoPhaseFluid[i].nuModel1().viscosityProperties().lookup("k3")
		)
	);

	k0Fluid2.set
	(
		i,
		new dimensionedScalar
		(
			"k0_2",
			dimensionSet(1, 1, -3, -1, 0, 0, 0),
			twoPhaseFluid[i].nuModel2().viscosityProperties().lookup("k0")
		)
	);

	k1Fluid2.set
	(
		i,
		new dimensionedScalar
		(
			"k1_2",
			dimensionSet(1, 1, -3, -2, 0, 0, 0),
			twoPhaseFluid[i].nuModel2().viscosityProperties().lookup("k1")
		)
	);

	k2Fluid2.set
	(
		i,
		new dimensionedScalar
		(
			"k2_2",
			dimensionSet(1, 1, -3, -3, 0, 0, 0),
			twoPhaseFluid[i].nuModel2().viscosityProperties().lookup("k2")
		)
	);

	k3Fluid2.set
	(
		i,
		new dimensionedScalar
		(
			"k3_2",
			dimensionSet(1, 1, -3, -4, 0, 0, 0),
			twoPhaseFluid[i].nuModel2().viscosityProperties().lookup("k3")
		)
	);

	Info<< "    Calculating phase 1 thermal conductivity kFluid1\n" << endl;
    kFluid1.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "k1",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            k0Fluid1[i] + k1Fluid1[i]*TFluid[i] + k2Fluid1[i]*TFluid[i]*TFluid[i] + k3Fluid1[i]*TFluid[i]*TFluid[i]*TFluid[i]
        )
    );

	Info<< "    Calculating phase 2 thermal conductivity kFluid2\n" << endl;
    kFluid2.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "k2",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            k0Fluid2[i] + k1Fluid2[i]*TFluid[i] + k2Fluid2[i]*TFluid[i]*TFluid[i] + k3Fluid2[i]*TFluid[i]*TFluid[i]*TFluid[i]
        )
    );

	Info<< "    Calculating average fluid thermal conductivity kFluid\n" << endl;
    kFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "k",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            alpha1l[i]*kFluid1[i] + alpha2l[i]*kFluid2[i]
        )
    );

	Info<< "    Adding to density rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            alpha1l[i]*rhoFluid1[i] + alpha2l[i]*rhoFluid2[i]
        )
    );

	Info<< "    Adding to specific heat cpFluid\n" << endl;
    cpFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "cp",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            alpha1l[i]*cpFluid1[i] + alpha2l[i]*cpFluid2[i]
        )
    );

	// Construct mass flux
	rhoPhiFluid.set
	(
		i,
		new surfaceScalarField
		(
			IOobject
			(
				"rhoPhi",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			fvc::interpolate(rhoFluid[i])*phiFluid[i]
		)
	);

	// Need to store rho*cp for ddt(rho*cp, T)
	// Initialization doesn't matter because rhoCp is reset after the gamma solution before use
	Info << "	Calculating rhoCp\n" << endl;
	rhoCpFluid.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"rhoCp",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			alpha1l[i]*rhoFluid1[i]*cpFluid1[i] 
			+ (scalar(1) - alpha1l[i])*rhoFluid2[i]*cpFluid2[i],
        	alpha1l[i].boundaryField().types()
		)
	);

	// Energy flux
	// Again, initialization doesn't matter as it is reset after gamma calculation
	Info << "	Calculating rhoPhiCp\n" << endl;
	rhoPhiCpFluid.set
	(
		i,
		new surfaceScalarField
		(
			IOobject
			(
				"rhoPhiCp",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			rhoPhiFluid[i]*cpFluid1[i]
		)
	);

	// Create incompressible turbulence model for each fluid region
    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
		autoPtr<incompressible::turbulenceModel>
		(
		    incompressible::turbulenceModel::New
		    (
		        UFluid[i],
		        phiFluid[i],
		        twoPhaseFluid[i]
			)
		).ptr()
    );

	Info<< "    Adding to gFluid\n" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, Zero) // uses name
        )
    );

    dimensionedScalar ghRef
    (
        mag(gFluid[i].value()) > SMALL
      ? gFluid[i]
          & (cmptMag(gFluid[i].value())/mag(gFluid[i].value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", gFluid[i].dimensions()*dimLength, 0)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (gFluid[i] & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (gFluid[i] & fluidRegions[i].Cf()) - ghRef
        )
    );

	Info<< "    Adding to p_rghFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

	Info<< "    Adding to pFluid\n" << endl;
    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
	Info<< "    Enforcing consistency between p and p_rgh\n" << endl;
    p_rghFluid[i] = pFluid[i] - rhoFluid[i]*ghFluid[i];

	Info<< "    Checking if flux required\n" << endl;
    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    Info<< "    Adding to radiationFluid\n" << endl;
	radiationFluid.set
    (
        i,
		autoPtr<radiation::radiationModel>
		(
		    radiation::radiationModel::New
		    (
		        TFluid[i]
			)
		).ptr()
    );

    const dictionary& pimpleDict = fluidRegions[i].solutionDict().subDict("PIMPLE");
    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);
	
	// Added to replace createDymControls.H
	pimpleDict.readIfPresent("correctPhi", correctPhiList[i]);
	pimpleDict.readIfPresent("checkMeshCourantNo", checkMeshCourantNoList[i]);
	pimpleDict.readIfPresent("moveMeshOuterCorrectors", moveMeshOuterCorrectorsList[i]);

    setRefCell
    (
        p_rghFluid[i],
        fluidRegions[i].solutionDict().subDict("PIMPLE"),
        pRefCellFluid[i],
        pRefValueFluid[i]
    );

	transportPropertiesl.set
	(
		i,
		new IOdictionary
		(
			IOobject
			(
				"transportProperties",
				runTime.constant(),
				fluidRegions[i],
				IOobject::MUST_READ_IF_MODIFIED,
				IOobject::NO_WRITE
			)
		)
	);

	// Create fields for traverse velocity
	//Info<< "	Creating velocity vector for modified reference frame\n" << endl;
	//VtFluid.set
	//(
	//	i,
	//	new dimensionedVector
	//	(
	//		"Vt",
	//		dimVelocity,
	//		transportPropertiesl[i].lookup("Vt")
	//	)
	//);

	UtFluid.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"Ut",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			fluidRegions[i],
			Vt //VtFluid[i]
		)
	);

	// Celcius conversion
	//celsConvl.set
	//(
	//	i,
	//	new dimensionedScalar
	//	(
	//		"celsConv",
	//		dimTemperature,
	//		transportPropertiesl[i].lookup("celsConv")
	//	)
	//);

	// Viscous heating coefficient
	phiViscl.set
	(
		i,
		new dimensionedScalar
		(
			"phiVisc",
			dimless,
			transportPropertiesl[i].lookup("phiVisc")
		)
	);

	// Create slip and frictional heat generation terms
	#include "createFrictionalHeating.H"

	// Create stress and strain variables if desired
	#include "createStressStrain.H"

	// Create particle tracking phase if desired
	#include "createBeta.H"

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    turbulenceFluid[i].validate();

	// MULES compressed flux left out, don't think it is needed (see interfoam)

}
