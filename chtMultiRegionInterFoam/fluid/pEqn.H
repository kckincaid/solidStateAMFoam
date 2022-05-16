{
    //if (correctPhi)
    //{
        //rAU.ref() = 1.0/UEqn.A();
    //}
    //else
    //{
        volScalarField rAU("rAU", 1.0/UEqn.A());
    //}

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    //volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    //if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        //adjustPhi(phiHbyA, U, p);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            thermo.surfaceTensionForce() //changed from mixture to thermo
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    //constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);
    constrainPressure(p, U, phiHbyA, rAUf, MRF);

    //while (pimple.correctNonOrthogonal())
	for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        //fvScalarMatrix pEqn
		fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
            //fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
        //pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
		//pEqn.solve(mesh.solver("pFinal"));

        //if (pimple.finalNonOrthogonalIter())
		if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA - p_rghEqn.flux();
            //phi = phiHbyA - pEqn.flux();

            p_rgh.relax();
			//p.relax();

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            //U = HbyA + rAU*fvc::reconstruct((phig - pEqn.flux())/rAUf);
			//U = Hby - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "incompressibleContinuityErrors.H"

	scalar maxU = max(mag(U)).value();
	reduce(maxU, maxOp<scalar>());

	Info << "Max velocity magnitude: " << maxU << " m/s" << endl;

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    //if (!correctPhi)
    //{
    //    rAU.clear();
    //}
}
