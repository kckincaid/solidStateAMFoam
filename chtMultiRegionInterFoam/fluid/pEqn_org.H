// didn't change from working chtMultiIco solver pEqn
// differences from interfoam: creation/use of rAUf face flux, solve for p_rgh
// unsure what effects these have

volScalarField rAU("rAU", 1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi))
    )
);

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
		fvScalarMatrix pEqn
		(
			fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
		);

		pEqn.setReference(pRefCell, pRefValue);

		pEqn.solve(mesh.solver("pFinal"));

        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    // Thermodynamic density update
    //thermo.correctRho(psi*p - psip0);
}

// Update continuity errors
#include "incompressibleContinuityErrors.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// Correct Uf if the mesh is moving
//fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
