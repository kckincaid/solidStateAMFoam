// Import appropriate mesh from region list
fvMesh& mesh = solidRegions[i];

// Set basic solid field values
volScalarField& temp = TSolid[i];
volScalarField& k = kSolid[i];
dimensionedScalar& rho = rhoSolid[i];
dimensionedScalar& cp = cpSolid[i];

// Create velocity and flux fields
volVectorField& U = USolid[i];
surfaceScalarField& phi = phiSolid[i];

// Variable thermal conductivity values
dimensionedScalar& k0s = k0Solid[i];
dimensionedScalar& k1s = k1Solid[i];
dimensionedScalar& k2s = k2Solid[i];
dimensionedScalar& k3s = k3Solid[i];

// Limiters for thermal conductivity
dimensionedScalar Tmax("Tmax", dimTemperature, solidProps[i].lookupOrDefault("Tmax", 1e+6));
dimensionedScalar Tmin("Tmin", dimTemperature, solidProps[i].lookupOrDefault("Tmin", 0));

// fvOptions for volumetric heat generation
fv::options& fvOptions = solidHeatSources[i];

// Update phi
phi = linearInterpolate(U) & mesh.Sf();
