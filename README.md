# solidStateAMFoam
Solid-state additive manufacturing within the OpenFOAM v1806 open-source CFD framework.

## Overview
This code was used to produce my PhD dissertation on solid-state additive manufacturing; specifically, simulating a novel solid-state process, additive friction stir deposition (AFSD). It was the first finite volume method model employed for this process, as only smoothed particle hydrodynamics methods had been used previously. The core of the code is formed by the solver and mixture model, with additional classes for the constitutive model, boundary conditions, and radiation model improvements. These components are described in brief below.

## Solver: chtMultiRegionInterFoam
The solver was based on the existing solver *chtMultiRegionFoam*, which is a coupled solver for heat transfer between an arbitrary number of solid and compressible fluid regions. Several steps were taken to modify this solver to represent the desired model physics, the first of which was to move from a fully compressible fluid description to a basic incompressible description, saving significant computational effort and being a quite acceptable approximation for this work. The second was to incorporate a mixture model to accurately represent both the build material and surrounding air; code from the two phase solver interFoam (using the Volume of Fluid method) was added to allow for this. Finally, dynamic mesh motion was incorporated to move the tool over the substrate mesh, if desired.

## Mixture class: immiscibleIncompressibleTwoPhaseThermalMixture
