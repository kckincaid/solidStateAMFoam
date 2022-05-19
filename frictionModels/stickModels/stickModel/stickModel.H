/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Namespace
    Foam::stickModels

Description
    A namespace for various stick models for solid-state material
	processing applications.

Class
    Foam::stickModel

Description
    An abstract base class for incompressible stickModels.

    The strain rate is defined by:

        mag(symm(grad(U)))


SourceFiles
    stickModel.C
    stickModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef stickModel_H
#define stickModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class stickModel Declaration
\*---------------------------------------------------------------------------*/

class stickModel
{

protected:

    // Protected data

        word name_;
        dictionary stickProperties_;

        const volVectorField& U_;
        const volScalarField& alpha_;
		const volScalarField& p_;


    // Private Member Functions

        //- No copy construct
        stickModel(const stickModel&) = delete;

        //- No copy assignment
        void operator=(const stickModel&) = delete;


public:

    //- Runtime type information
    TypeName("stickModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            stickModel,
            dictionary,
            (
                const word& name,
                const dictionary& stickProperties,
                const volVectorField& U,
                const volScalarField& alpha,
				const volScalarField& p
            ),
            (name, stickProperties, U, alpha, p)
        );


    // Selectors

        //- Return a reference to the selected stick model
        static autoPtr<stickModel> New
        (
            const word& name,
            const dictionary& stickProperties,
            const volVectorField& U,
            const volScalarField& alpha,
			const volScalarField& p
        );


    // Constructors

        //- Construct from components
        stickModel
        (
            const word& name,
            const dictionary& stickProperties,
            const volVectorField& U,
            const volScalarField& alpha,
			const volScalarField& p
        );


    //- Destructor
    virtual ~stickModel()
    {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& stickProperties() const
        {
            return stickProperties_;
        }

        //- Return the strain rate
        tmp<volScalarField> strainRate() const;

        //- Return the internal viscous heating
        virtual tmp<volScalarField> qvisc() const = 0;

		//- Return the cell-centered frictional heating
		virtual tmp<surfaceScalarField> qfric() const = 0;

        //- Correct the heating quantities
        virtual void correct() = 0;

        //- Read frictionProperties dictionary
        virtual bool read(const dictionary& stickProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
