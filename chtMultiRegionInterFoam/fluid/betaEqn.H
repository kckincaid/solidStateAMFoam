{
	// Solve for beta, the "particle tracking" field
	tmp<fvScalarMatrix> tBetaEqn
	(
		fvm::ddt(rho, beta)
	  + fvm::div(rhoPhi, beta)
	  - fvm::laplacian(rho*gammaBeta, beta)
	);

	fvScalarMatrix& betaEqn = tBetaEqn.ref();

	betaEqn.solve(mesh.solver(beta.select(finalIter)));
}
