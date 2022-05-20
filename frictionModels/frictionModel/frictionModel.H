/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::frictionModel

Description
    Base-class for all friction models.

SourceFiles
    frictionModel.C

\*---------------------------------------------------------------------------*/

#ifndef frictionModel_H
#define frictionModel_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class frictionModel Declaration
\*---------------------------------------------------------------------------*/

class frictionModel
{
    // Private Member Functions

        //- No copy construct
        frictionModel(const frictionModel&) = delete;

        //- No copy assignment
        void operator=(const frictionModel&) = delete;


public:

    //- Runtime type information
    TypeName("frictionModel");


    // Constructors

        //- Construct from components
        frictionModel
        ();


    //- Destructor
    virtual ~frictionModel();


    // Member Functions

        //- Return the viscous heating
        virtual tmp<volScalarField> qvisc() const = 0;

        //- Return the cell-centered frictional heating
        virtual tmp<volScalarField> qfric() const = 0;

        //- Correct the frictional and viscous heating
        virtual void correct() = 0;

        //- Read frictionProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
