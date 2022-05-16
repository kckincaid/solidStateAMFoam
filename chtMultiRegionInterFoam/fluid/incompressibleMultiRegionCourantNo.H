    scalar CoNum = -GREAT;

    forAll(fluidRegions, regioni)
    {
        CoNum = max
        (
            incompressibleCourantNo
            (
                fluidRegions[regioni],
                runTime,
                phiFluid[regioni]
            ),
            CoNum
        );
    }
/*
    forAll(porousFluidRegions, porousi)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                porousFluidRegions[porousi],
                runTime,
                rhoPorous[porousi],
                phiPorous[porousi]
            ),
            CoNum
        );
    }
*/
