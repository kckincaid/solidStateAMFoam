    // Solve the Momentum equation
	// Copied mostly from interfoam UEqn.H file

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U) 
	  + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
			//-fvc::grad(p) // from pimple/icofoam
            fvc::reconstruct
            (
                (
                    thermo.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
				  //- fvc::snGrad(p)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
