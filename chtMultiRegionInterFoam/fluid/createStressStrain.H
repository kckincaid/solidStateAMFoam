// Check if stress / strain calculations are desired
transportPropertiesl[i].readIfPresent("trackStressStrain", trackStressStrainl[i]);

// Read in variables if strain tracking desired
if (trackStressStrainl[i])
{
    Info<< "	Creating stress and strain for region " << fluidRegions[i].name() << nl << endl;

    // Zener-Hollomon parameter
    ZFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Z",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
			dimensionedScalar(dimless/dimTime, Zero)
        )
    );

    // Flow stress
    sigmafFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "sigmaf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
			dimensionedScalar(dimPressure, Zero)
        )
    );

    // Equivalent plastic strain
    epsilonEqFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "epsilonEq",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
			dimensionedScalar(dimless, Zero)
        )
    );
    // Diffusion coefficient for strain transport eqn
    gammaEpsilonFluid.set
    (
        i,
        new dimensionedScalar
        (
            "gammaEpsilon",
            dimLength*dimLength/dimTime,
            transportPropertiesl[i].lookup("gammaEpsilon")
        )
    );
    // von Mises stress
    sigmavmFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "sigmavm",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimPressure, Zero)
        )
    );
}

// Create dummy variables if not needed
else
{
    // Zener-Hollomon parameter
    ZFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Z",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i],
			dimensionedScalar(dimless/dimTime, Zero)
        )
    );

    // Flow stress
    sigmafFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "sigmaf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i],
			dimensionedScalar(dimPressure, Zero)
        )
    );

    // Equivalent plastic strain
    epsilonEqFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "epsilonEq",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
    // Diffusion coefficient for strain transport eqn
    gammaEpsilonFluid.set
    (
        i,
        new dimensionedScalar
        (
            "gammaEpsilon",
            dimLength*dimLength/dimTime,
            Zero
        )
    );
    // von Mises stress
    sigmavmFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "sigmavm",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimPressure, Zero)
        )
    );
}
