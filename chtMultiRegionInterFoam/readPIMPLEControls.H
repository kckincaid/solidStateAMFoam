    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    fvSolution solutionDict(runTime);

    const dictionary& pimpleDict = solutionDict.subDict("PIMPLE");

    const int nOuterCorr =
        pimpleDict.lookupOrDefault<int>("nOuterCorrectors", 1);
