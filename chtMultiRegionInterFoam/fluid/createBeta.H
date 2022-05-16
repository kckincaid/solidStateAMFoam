// Check if particle tracking is desired
transportPropertiesl[i].readIfPresent("particleTracking", particleTrackingl[i]);

// Create beta field and read in diffusion coefficient if desired
if (particleTrackingl[i])
{
    Info<< "Reading tracking phase properties for region " << fluidRegions[i].name() << nl << endl;

    betaFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "beta",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );
    // Phase diffusion coefficient
    gammaBetaFluid.set
    (
        i,
        new dimensionedScalar
        (
            "gammaBeta",
            dimLength*dimLength/dimTime,
            transportPropertiesl[i].lookup("gammaBeta")
        )
    );
}

// Create dummy fields if not desired
else
{
    betaFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "beta",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
    // Phase diffusion coefficient
    gammaBetaFluid.set
    (
        i,
        new dimensionedScalar
        (
            "gammaBeta",
            dimLength*dimLength/dimTime,
            Zero
        )
    );
}
