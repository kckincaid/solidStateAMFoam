// Create new dictionary in constant for fields used throughout solver
IOdictionary caseProperties
(
	IOobject
	(
		"caseProperties",
		runTime.constant(),
		runTime.db(),
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

// Read in moving reference frame velocity vector from caseProperties
dimensionedVector Vt("Vt", dimVelocity, caseProperties.lookup("Vt"));

// Read in interpolation method for calculating
const word interpolationType = caseProperties.lookupOrDefault<word>("interpolationType", "arithmetic");

#include "createFluidFields.H"
#include "createSolidFields.H"
