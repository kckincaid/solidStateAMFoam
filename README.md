# solidStateAMFoam
Simulating solid-state additive manufacturing within the OpenFOAM (v1806) open-source CFD framework.

## Overview
This code was used to produce my PhD dissertation on solid-state additive manufacturing; specifically, simulating a novel solid-state process, additive friction stir deposition (AFSD). It was the first finite volume method model employed for this process, as only smoothed particle hydrodynamics methods had been used previously. The core of the code is formed by the solver and mixture model, with additional classes for the constitutive and friction models, boundary conditions, and radiation model improvements. These components are described in brief below.

## Solver: chtMultiRegionInterFoam
The solver was based on the existing solver *chtMultiRegionFoam*, which is a coupled solver for heat transfer between an arbitrary number of solid and compressible fluid regions. Several steps were taken to modify this solver to represent the desired model physics, the first of which was to move from a fully compressible fluid description to a basic incompressible description, saving significant computational effort and being a quite acceptable approximation for this work. The second was to incorporate a mixture model to accurately represent both the build material and surrounding air; code from the two phase solver interFoam (using the Volume of Fluid method) was added to allow for this. Finally, dynamic mesh motion was incorporated to move the tool over the substrate mesh, if desired.

## Mixture class: immiscibleIncompressibleTwoPhaseThermalMixture
The mixture class is a critical accessory to the solver, calculating the thermophysical properties of the fluid throughout the domain as time progresses. The density, specific heat, and thermal conductivity are directly updated from this class, and the viscosity update is called through it. It was based on the immiscibleIncompressibleTwoPhaseMixture class, with additions for the thermal properties. Currently it supports constant density and specific heat, as well as either a constant thermal conductivity or a cubic fit temperature-dependent conductivity with user-supplied coefficients.

## Viscoplastic models
The viscoplastic constitutive models are responsible for relating the effective viscosity of the build material to its state. Specifically, three models are implemented in which the viscosity depends on the strain rate and temperature of the material. The first of these is the Sellars-Tegart (a.k.a. Sheppard-Wright) model, which has been extensively used in other hot working applications (e.g. friction stir welding) and was shown to provide excellent results in this application as well. The second is based off of the Fraser-Kiss-St. Georges (FKS) model, which typically includes a hardening effect based on the plastic strain. The implementation herein assumes a fully saturated plastic strain (epsilon > 0.1), which is reasonable for SS AM applications. The final implemented model is Norton-Hoff with constant coefficients, which is largely untested at this point and is expected to provide poor results.

## Friction model
This class calculates the viscous and frictional heating terms present in solid-state processes. Two model types are presently implemented, differing in the slip condition prescribed at the interface. Both employ a slip/stick assumption, but one assumes a constant slip fraction while the other prescribes a slip value dependent on the distance from the tool center. The viscous heating (actually a result of plastic dissipation) is calculted from the strain rate and effective viscosity of the material, and generally remains the same between models.

## Boundary conditions
A suite of boundary conditions has been implemented to reduce dependence on GroovyBC, a library included in the swak4FOAM add-on package. These include conditions for both heat transfer (temperature) as well as velocity for a number of applications.

## Radiation
The solver incorporates radiation for solving heat transfer in extreme conditions. A two phase absorption / emission model was created to work with the VOF solver to allow different radiative properties for the two phases.
