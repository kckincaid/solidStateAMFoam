/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::incompressibleFrictionModel

Description
    A simple friction model based on stickModel.

SourceFiles
    incompressibleFrictionModel.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleFrictionModel_H
#define incompressibleFrictionModel_H

#include "frictionModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class stickModel;

/*---------------------------------------------------------------------------*\
                Class incompressibleFrictionModel Declaration
\*---------------------------------------------------------------------------*/

class incompressibleFrictionModel
:
    public IOdictionary,
    public frictionModel
{
    // Private Data

        autoPtr<stickModel> stickModelPtr_;


    // Private Member Functions

        //- No copy construct
        incompressibleFrictionModel(const incompressibleFrictionModel&) = delete;

        //- No copy assignment
        void operator=(const incompressibleFrictionModel&) = delete;


public:

    //- Runtime type information
    TypeName("incompressibleFrictionModel");


    // Constructors

        //- Construct from components
        incompressibleFrictionModel
        (
            const volVectorField& U,
			const volVectorField& Ut
        );


    //- Destructor
    virtual ~incompressibleFrictionModel();


    // Member Functions

        //- Return the internal viscous heating
        virtual tmp<volScalarField> qvisc() const;

        //- Return the cell-center frictional heating
        virtual tmp<volScalarField> qfric() const;

        //- Correct the heating quantities
        virtual void correct();

        //- Read frictionProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
