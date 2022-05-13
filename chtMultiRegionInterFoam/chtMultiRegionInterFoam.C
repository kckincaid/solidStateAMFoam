/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtMultiRegionInterFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for incompressible, turbulent fluid flow with frictional
	heat generation, and solid heat conduction with conjugate heat transfer 
	between solid and fluid regions.

    It handles secondary fluid or solid circuits which can be coupled
    thermally with the main fluid region. i.e radiators, etc.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"				// from pimplefoam - added for dym
//#include "turbulentFluidThermoModel.H"
//#include "rhoReactionThermo.H"
//#include "CombustionModel.H"
#include "singlePhaseTransportModel.H" 	// from pimplefoam
#include "turbulentTransportModel.H" 	// from pimplefoam
#include "fixedGradientFvPatchFields.H" // from pimplefoam
#include "pimpleControl.H"				// from pimplefoam - added for dym
#include "CorrectPhi.H"					// from pimplefoam
#include "CMULES.H"						// from interfoam
#include "EulerDdtScheme.H"				// from interfoam
#include "localEulerDdtScheme.H"		// from interfoam
#include "CrankNicolsonDdtScheme.H"		// from interfoam
#include "subCycle.H"					// from interfoam
#include "immiscibleIncompressibleTwoPhaseThermalMixture.H" // new class
#include "regionProperties.H"
#include "incompressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //#define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
	// PimpleFoam has "createDymanicFvMesh.H" - look at it and modify the
	// createMesh.H files herein to match them.
    #include "createMeshes.H"
	//#include "createDyMControls.H" - moved into createFluidFields and setRegionFluidFields
    #include "createFields.H"
	#include "createAlphaFluxes.H" 			// from interfoam
    #include "initContinuityErrs.H"
    #include "createTimeControls.H" // pre included in createDymControls
    #include "readSolidTimeControls.H"
    #include "incompressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    while (runTime.run())
    {
		//#include "readDyMControls.H"
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        #include "incompressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldFluidFields.H"
            }
        }

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)
        {
            const bool finalIter = (oCorr == nOuterCorr-1);

            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "solveSolid.H"
            }

            // Additional loops for energy solution only
            if (!oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimpleDict, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                       #include "setRegionFluidFields.H"
                       #include "readFluidMultiRegionPIMPLEControls.H"
                       frozenFlow = true;
                       #include "solveFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "setRegionSolidFields.H"
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "solveSolid.H"
                    }
                }
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
