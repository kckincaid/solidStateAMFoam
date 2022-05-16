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

	// Mixture thermal properties
	volScalarField& rho = rhoFluid[i];
	volScalarField& kf = kFluid[i];
	volScalarField& cpf = cpFluid[i];

	// Individual phase properties
	dimensionedScalar rho1 = thermo.rho1();
	dimensionedScalar rho2 = thermo.rho2();
	dimensionedScalar cp1	= thermo.cp1();
	dimensionedScalar cp2  = thermo.cp2();

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
