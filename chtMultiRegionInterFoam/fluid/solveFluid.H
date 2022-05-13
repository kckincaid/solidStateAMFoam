if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (pimple.firstIter() || moveMeshOuterCorrectors)
{
	mesh.update();

	if (mesh.changing())
    {
		// Do not apply previous time-step mesh compression flux
        // if the mesh topology changed (from interfoam)
        if (mesh.topoChanging())
        {
            talphaPhi1Corr0.clear();
        }

        MRF.update();

        if (correctPhi)
        {
            // Calculate absolute flux
            // from the mapped surface velocity
            phi = mesh.Sf() & Uf();

            #include "correctPhi.H"

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);

			thermo.correct(); // from interfoam
        }

        if (checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }
    }
}

if (frozenFlow)
{
    #include "TEqn.H"
}
else
{	
	// Solve for phase fraction alpha
	#include "alphaControls.H"
	#include "alphaEqnSubCycle.H"

	thermo.correct();

    #include "UEqn.H"
    #include "TEqn.H"

    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"
    }

    turbulence.correct();

}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}

// Solve for strain and tracked particles after velocity field is finalized
if (trackStressStrain)
{
    #include "calculateStrain.H"
}

if (particleTracking)
{
    #include "betaEqn.H"
}

