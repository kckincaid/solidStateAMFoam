    dynamicFvMesh& mesh = fluidRegions[i];

	// Include this to define LTS variable after mesh defined
	#include "createRDeltaT.H"

	immiscibleIncompressibleTwoPhaseThermalMixture& thermo = twoPhaseFluid[i];

	// Create basic fluid property fields
	volScalarField& p = pFluid[i];
	volScalarField& p_rgh = p_rghFluid[i];
	volScalarField& temp = TFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

	// Create gravity fields for pressure calcs
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

	// Viscosity and viscous heating terms
	dimensionedScalar& phiVisc = phiViscl[i];
	volScalarField& muEff = muEffl[i];
	volScalarField& epsilonDotEq = epsilonDotEql[i];
	volScalarField& qvisc = qviscl[i];

    incompressible::turbulenceModel& turbulence = turbulenceFluid[i];

	volScalarField& rho = rhoFluid[i];
	volScalarField& kf = kFluid[i];

	dimensionedScalar& rho1 = rhoFluid1[i];
	dimensionedScalar& rho2 = rhoFluid2[i];
	dimensionedScalar& cp1	= cpFluid1[i];
	dimensionedScalar& cp2 = cpFluid2[i];

	//volScalarField& k1f = kFluid1[i];
	//volScalarField& k2f = kFluid2[i];

	// Temperature-dependent thermal conductivity constants
	//dimensionedScalar& k0_1 = k0Fluid1[i];
	//dimensionedScalar& k1_1 = k1Fluid1[i];
	//dimensionedScalar& k2_1 = k2Fluid1[i];
	//dimensionedScalar& k3_1 = k3Fluid1[i];
	//dimensionedScalar& k0_2 = k0Fluid2[i];
	//dimensionedScalar& k1_2 = k1Fluid2[i];
	//dimensionedScalar& k2_2 = k2Fluid2[i];
	//dimensionedScalar& k3_2 = k3Fluid2[i];

	// Temperature-dependent max shear stress constants
	dimensionedScalar& TC0 = TC0Fluid[i];
	dimensionedScalar& TC1 = TC1Fluid[i];
	dimensionedScalar& TC2 = TC2Fluid[i];
	dimensionedScalar& TC3 = TC3Fluid[i];
	dimensionedScalar& TC4 = TC4Fluid[i];
	dimensionedScalar& Tmin = TminFluid[i];
	dimensionedScalar& Tmax = TmaxFluid[i];
	dimensionedScalar& taumin = tauminFluid[i];
	dimensionedScalar& taumax = taumaxFluid[i];

	// Stress and strain fields and variables
	bool trackStressStrain = trackStressStrainl[i];
	volScalarField& Z = ZFluid[i];
	volScalarField& sigmaf = sigmafFluid[i];
	volScalarField& epsilonEq = epsilonEqFluid[i];
	dimensionedScalar& gammaEpsilon = gammaEpsilonFluid[i];
	volScalarField& sigmavm = sigmavmFluid[i];

	// Particle tracking / 2nd phase
	bool particleTracking = particleTrackingl[i];
	volScalarField& beta = betaFluid[i];
	dimensionedScalar& gammaBeta = gammaBetaFluid[i];
	
	//dimensionedScalar& celsConv = celsConvl[i];

    volScalarField& alpha1(thermo.alpha1());
	volScalarField& alpha2(thermo.alpha2());

	surfaceScalarField& rhoPhi = rhoPhiFluid[i];
	surfaceScalarField& rhoPhiCp = rhoPhiCpFluid[i];
	volScalarField& rhoCp = rhoCpFluid[i];

    radiation::radiationModel& rad = radiationFluid[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    bool frozenFlow = frozenFlowFluid[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

	// Alpha solver variables
	bool& alphaRestart = alphaRestartList[i];
	surfaceScalarField& alphaPhi10 = alphaPhi10List[i];
	tmp<surfaceScalarField> talphaPhi1Corr0;

	#include "createUfIfPresent.H"

	// Frictional heat generation terms
	const wordList& fricPatchList = transportPropertiesl[i].lookup("fricPatch");
	const word& slipType = slipTypel[i];
	const dimensionedScalar& etaf = etafl[i];
	const dimensionedScalar& omega = omegal[i];
	const dimensionedScalar& omega0 = omega0l[i];
	const dimensionedScalar& delta0 = delta0l[i];
	const dimensionedScalar& Rs = Rsl[i];
	const dimensionedScalar& mu0 = mu0l[i];
	const dimensionedScalar& lambda = lambdal[i];
	volScalarField& qfric = qfricl[i];

	// Patch center locating fields - can remove when patch center calculator works
	const dimensionedVector& startPos = transportPropertiesl[i].lookup("startPos");
	const dimensionedVector& patchVel = transportPropertiesl[i].lookup("patchVel");

	// Dynamic Mesh controls
	bool& correctPhi = correctPhiList[i];
	bool& checkMeshCourantNo = checkMeshCourantNoList[i];
	bool& moveMeshOuterCorrectors = moveMeshOuterCorrectorsList[i];

	pimpleControl pimple(mesh);

	#include "initCorrectPhi.H"
