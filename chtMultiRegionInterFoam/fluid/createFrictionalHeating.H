// Look up type of condition to apply at slip faces
transportPropertiesl[i].readIfPresent("slipType", slipTypel[i]);

// If constant slip selected, read in the following:
//      omega       rotational speed
//      delta0      slip fraction
//      mu0         friction coefficient
//      Rs          tool shoulder radius
//      etaf        mechanical efficiency
// and set the following variables to dummy values (do not read):
//      omega0
//      lambda
if (slipTypel[i] == "constant")
{
    // Set patch rotational speed
    omegal.set
    (
        i,
        new dimensionedScalar
        (
            "omega",
            dimless/dimTime,
            transportPropertiesl[i].lookup("omega")
        )
    );

    // Set slip fraction on patch face
    delta0l.set
    (
        i,
        new dimensionedScalar
        (
            "delta0",
            dimless,
            transportPropertiesl[i].lookup("delta0")
        )
    );

    // Set friction coefficient on patch face
    mu0l.set
    (
        i,
        new dimensionedScalar
        (
            "mu0",
            dimless,
            transportPropertiesl[i].lookup("mu0")
        )
    );

    // Set patch radius
    Rsl.set
    (
        i,
        new dimensionedScalar
        (
            "Rs",
            dimLength,
            transportPropertiesl[i].lookup("Rs")
        )
    );

    // Set mechanical efficiency of plastic work
    etafl.set
    (
        i,
        new dimensionedScalar
        (
            "etaf",
            dimless,
            transportPropertiesl[i].lookup("etaf")
        )
    );

    // Dummy variable: exponential slip coefficient
    omega0l.set
    (
        i,
        new dimensionedScalar
        (
            "omega0",
            dimless/dimTime,
            Zero
        )
    );

    // Dummy variable: exponential friction coefficient
    lambdal.set
    (
        i,
        new dimensionedScalar
        (
            "lambda",
            dimless/dimVelocity,
            Zero
        )
    );
}

// If slip type set to "exponential," read in all variables
else if (slipTypel[i] == "exponential")
{
    // Set patch rotational speed
    omegal.set
    (
        i,
        new dimensionedScalar
        (
            "omega",
            dimless/dimTime,
            transportPropertiesl[i].lookup("omega")
        )
    );

    // Set slip fraction on patch face
    delta0l.set
    (
        i,
        new dimensionedScalar
        (
            "delta0",
            dimless,
            transportPropertiesl[i].lookup("delta0")
        )
    );

    // Set friction coefficient on patch face
    mu0l.set
    (
        i,
        new dimensionedScalar
        (
            "mu0",
            dimless,
            transportPropertiesl[i].lookup("mu0")
        )
    );

    // Set patch radius
    Rsl.set
    (
        i,
        new dimensionedScalar
        (
            "Rs",
            dimLength,
            transportPropertiesl[i].lookup("Rs")
        )
    );

    // Set mechanical efficiency of plastic work
    etafl.set
    (
        i,
        new dimensionedScalar
        (
            "etaf",
            dimless,
            transportPropertiesl[i].lookup("etaf")
        )
    );

    // Set exponential slip coefficient
    omega0l.set
    (
        i,
        new dimensionedScalar
        (
            "omega0",
            dimless/dimTime,
            transportPropertiesl[i].lookup("omega0")
        )
    );

    // Dummy variable: exponential friction coefficient
    lambdal.set
    (
        i,
        new dimensionedScalar
        (
            "lambda",
            dimless/dimVelocity,
            transportPropertiesl[i].lookup("lambda")
        )
    );
}

// If constant slip selected, read in the following:
//      omega       rotational speed
//      mu0         friction coefficient
//      Rs          tool shoulder radius
//      etaf        mechanical efficiency
// and set the following variables to dummy values (do not read):
//		delta0
//      omega0
//      lambda
else if (slipTypel[i] == "friction")
{
    // Set patch rotational speed
    omegal.set
    (
        i,
        new dimensionedScalar
        (
            "omega",
            dimless/dimTime,
            transportPropertiesl[i].lookup("omega")
        )
    );

    // Set friction coefficient on patch face
    mu0l.set
    (
        i,
        new dimensionedScalar
        (
            "mu0",
            dimless,
            transportPropertiesl[i].lookup("mu0")
        )
    );

    // Set patch radius
    Rsl.set
    (
        i,
        new dimensionedScalar
        (
            "Rs",
            dimLength,
            transportPropertiesl[i].lookup("Rs")
        )
    );

    // Set mechanical efficiency of plastic work
    etafl.set
    (
        i,
        new dimensionedScalar
        (
            "etaf",
            dimless,
            transportPropertiesl[i].lookup("etaf")
        )
    );


    // Dummy variable: slip fraction
    delta0l.set
    (
        i,
        new dimensionedScalar
        (
            "delta0",
            dimless,
            Zero
        )
    );

    // Dummy variable: exponential slip coefficient
    omega0l.set
    (
        i,
        new dimensionedScalar
        (
            "omega0",
            dimless/dimTime,
            Zero
        )
    );

    // Dummy variable: exponential friction coefficient
    lambdal.set
    (
        i,
        new dimensionedScalar
        (
            "lambda",
            dimless/dimVelocity,
            Zero
        )
    );
}

// Throw error if slip type not recognized
else
{
    FatalErrorIn
        (
            "Slip function not recognized"
        )   << "Choose constant, exponential, or friction slip"
            << exit(FatalError);
}

// Temperature-dependent max shear stress
TC0Fluid.set
(
    i,
    new dimensionedScalar
    (
        "TC0",
        dimPressure,
        transportPropertiesl[i].lookup("TC0")
    )
);

TC1Fluid.set
(
    i,
    new dimensionedScalar
    (
        "TC1",
        dimPressure,
        transportPropertiesl[i].lookup("TC1")
    )
);

TC2Fluid.set
(
    i,
    new dimensionedScalar
    (
        "TC2",
        dimPressure,
        transportPropertiesl[i].lookup("TC2")
    )
);

TC3Fluid.set
(
    i,
    new dimensionedScalar
    (
        "TC3",
        dimPressure,
        transportPropertiesl[i].lookup("TC3")
    )
);

TC4Fluid.set
(
    i,
    new dimensionedScalar
    (
        "TC4",
        dimPressure,
        transportPropertiesl[i].lookup("TC4")
    )
);

TminFluid.set
(
    i,
    new dimensionedScalar
    (
        "Tmin",
        dimTemperature,
        transportPropertiesl[i].lookup("Tmin")
    )
);

TmaxFluid.set
(
    i,
    new dimensionedScalar
    (
        "Tmax",
        dimTemperature,
        transportPropertiesl[i].lookup("Tmax")
    )
);

tauminFluid.set
(
    i,
    new dimensionedScalar
    (
        "taumin",
        dimPressure,
        transportPropertiesl[i].lookup("taumin")
    )
);

taumaxFluid.set
(
    i,
    new dimensionedScalar
    (
        "taumax",
        dimPressure,
        transportPropertiesl[i].lookup("taumax")
    )
);

// Create viscous heating field
qfricl.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "qfric",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluidRegions[i]
    )
);
