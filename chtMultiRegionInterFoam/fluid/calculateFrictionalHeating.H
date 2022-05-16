	// Reset frictional heat term
	qfric = qfric*0.0;

	// Calculate alpha phase fraction on faces
	surfaceScalarField alphaf = fvc::interpolate(alpha1);

	// Evaluate delta, muf, and qfric for cells adjacent to frictional heating patches
	// Note: fricPatchList is user-defined list of patches where heating should be applied
	forAll(fricPatchList, patchi)
	{
		// Get current patch
		scalar patchNum = mesh.boundaryMesh().findPatchID(fricPatchList[patchi]);
		const polyPatch& currPatch = mesh.boundaryMesh()[patchNum];

		// Initialize center vector and area
		vector patchCenter = Zero;
		scalar patchArea = VSMALL;

		// Create pointers to face surface area and centers (used by both loops)
		const surfaceScalarField& magSf = mesh.magSf();
		const surfaceVectorField& Cf = mesh.Cf();

		// Calculate area-weighted patch center
		forAll(currPatch, facei)
		{
			// Add area-weighted face center and face area to totals
			patchCenter += Cf.boundaryField()[patchNum][facei]*magSf.boundaryField()[patchNum][facei];
			patchArea += magSf.boundaryField()[patchNum][facei];
		}

		reduce(patchCenter, sumOp<vector>());
		reduce(patchArea, sumOp<scalar>());

		// Calculate patch centroid
		patchCenter = patchCenter / patchArea;

		// Workaround for now: calculate patch center using old method
		//patchCenter = patchVel.value()*runTime.value() + startPos.value();

		Info << "Center for patch " << fricPatchList[patchi] << " found to be at " << patchCenter << endl;

		// Initialize average radius (just for sanity check)
		//scalar avgRad = 0;
		scalar patchSize = 0;

		// Loop over all faces on patch		
		forAll(currPatch, facei)
		{
			// Get labels for cells which own patch faces
			label faceCelli = currPatch.faceCells()[facei];

			// Get distance of cell center from center of rotation
			scalar dx = Cf.boundaryField()[patchNum][facei].component(0) - patchCenter.component(0);
			scalar dy = Cf.boundaryField()[patchNum][facei].component(1) - patchCenter.component(1);
			scalar dz = Cf.boundaryField()[patchNum][facei].component(2) - patchCenter.component(2);

			scalar rad = Foam::sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));

			//avgRad += rad;
			patchSize++;

			//dimensionedScalar rad = Foam::sqrt(pow(mesh.C()[faceCelli].component(0) - startX.value()
											  	//- runTime.value()*depVel.value(),2) 
											  //+pow(mesh.C()[faceCelli].component(1),2))*Rs/Rs.value();
			
			// Calculate angle of point w.r.t origin
			//dimensionedScalar theta = Foam::atan(mesh.C()[faceCelli].component(1)/mesh.C()[faceCelli].component(2));
			scalar theta = Foam::atan(Cf.boundaryField()[patchNum][facei].component(1) / Cf.boundaryField()[patchNum][facei].component(0));

			// Create dummy variables for delta (slip) and mu (friction coefficient)
			dimensionedScalar delta = 1;
			dimensionedScalar muf = 1;

			// If constant slip, set to specified value
			if (slipType == "constant")
			{
				delta = delta0;
				muf = mu0;
			}

			// If exponential slip, calculate slip and muf at cell centers
			if (slipType == "exponential")
			{
				delta -= exp(-(omega*rad)/(delta0*omega0*Rs.value()));
				muf = mu0*exp(-lambda.value()*delta*omega.value()*rad);
			}

			// If friction slip, calculate slip from material U and patch U
			if (slipType == "friction")
			{
				// Calculate patch velocity field
				dimensionedScalar Up = rad*omega;

				// Calculate local slip fraction
				delta = mag(U[i])/Up.value();

				// Set friction coefficient
				muf = mu0;
			}

			// Calculate temperature-dependent yield shear stress
			dimensionedScalar tau = TC0*0.0;

			if(temp[faceCelli] < Tmin.value())
			{
				tau = taumax;
			}
			else if(temp[faceCelli] > Tmax.value())
			{
				tau = taumin;
			}
			else
			{
				// Calculate temp-dependent max shear stress if temperature is within bounds
				// However correlation is degrees K to MPa, so multiply by 1e+6 to convert to Pa
				tau = (TC0 + TC1*temp[faceCelli] + TC2*pow(temp[faceCelli],2)
					+ TC3*pow(temp[faceCelli],3) + TC4*pow(temp[faceCelli],4))*1e+6;
			}

			// Generate frictional heating field
			qfric[faceCelli] = alpha1[faceCelli]*((scalar(1) - delta.value())*etaf.value()*tau.value()
								+ delta.value()*muf.value()*p_rgh[faceCelli])*(omega.value()*rad 
								- Vt.component(0).value()*Foam::sin(theta));
		}

		// Calculate average radius (for sanity check)
		//avgRad = avgRad / patchSize;
		//Info << "Patch size " << patchSize << endl;

		//Info << "Average radius for patch " << fricPatchList[patchi] << " found to be " << avgRad << endl;
	}

	qfric.correctBoundaryConditions();
