/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "functionObjectTemplate.H"
#define namespaceFoam  // Suppress <using namespace Foam;>
#include "fvCFD.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(calcMomentumFluxFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    calcMomentumFluxFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = c30f89b3a91cbed40703f29b7c17826632a6186c
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void calcMomentumFlux_c30f89b3a91cbed40703f29b7c17826632a6186c(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fvMesh&
Foam::calcMomentumFluxFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
calcMomentumFluxFunctionObject::
calcMomentumFluxFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
calcMomentumFluxFunctionObject::
~calcMomentumFluxFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
calcMomentumFluxFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read calcMomentumFlux");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcMomentumFluxFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute calcMomentumFlux");
    }

//{{{ begin code
    #line 119 "/home/enrique/Downloads/casos_1_2/caso_corregido_165/system/controlDict/functions/inletMomentumFlux"
// Get mesh and fields
            const fvMesh& mesh = refCast<const fvMesh>(obr_);
            const volVectorField& U = mesh.lookupObject<volVectorField>("U");
            const volScalarField& alpha = mesh.lookupObject<volScalarField>("alpha.fuel");
            const surfaceScalarField& phi = mesh.lookupObject<surfaceScalarField>("phi");
            const volScalarField& p = mesh.lookupObject<volScalarField>("p_rgh");
            
            // Density values
            scalar rhoFuel = 1000.0;
            scalar rhoAir = 20.0;
            
            // Wedge correction factor
            scalar wedgeFactor = 72.0;
            
            // Get inlet patch
            label patchI = mesh.boundaryMesh().findPatchID("inlet");
            const fvPatch& patch = mesh.boundary()[patchI];
            
            // Get face values
            const vectorField& Sf = patch.Sf();
            const scalarField& magSf = patch.magSf();
            const scalarField& phip = phi.boundaryField()[patchI];
            const vectorField& Up = U.boundaryField()[patchI];
            const scalarField& alphap = alpha.boundaryField()[patchI];
            const scalarField& pp = p.boundaryField()[patchI];
            
            // Calculate forces
            scalar massFlow = 0.0;
            vector momentumFlux = vector::zero;
            vector pressureForce = vector::zero;
            scalar inletArea = 0.0;
            
            forAll(patch, faceI)
            {
                scalar rho = alphap[faceI]*rhoFuel + (1.0 - alphap[faceI])*rhoAir;
                scalar mDot = rho * phip[faceI];
                massFlow += mDot;
                momentumFlux += mDot * Up[faceI];
                pressureForce += pp[faceI] * Sf[faceI];
                inletArea += magSf[faceI];
            }
            
            // Apply wedge correction
            scalar massFlowFull = massFlow * wedgeFactor;
            vector momentumFluxFull = momentumFlux * wedgeFactor;
            vector pressureForceFull = pressureForce * wedgeFactor;
            scalar inletAreaFull = inletArea * wedgeFactor;
            
            // Total force (momentum + pressure)
            vector totalForce = momentumFluxFull + pressureForceFull;
            
            // Output file path
            fileName outputFile = mesh.time().path()/"postProcessing"/"jetForces.dat";
            
            // Create postProcessing directory if needed
            mkDir(mesh.time().path()/"postProcessing");
            
            // Check if file exists using stat
            struct stat buffer;
            bool fileExists = (stat(outputFile.c_str(), &buffer) == 0);
            
            // Open file in append mode
            std::ofstream file(outputFile.c_str(), std::ios::app);
            
            if (!fileExists)
            {
                file << "# Jet Force Analysis - OpenFOAM Post-Processing\n";
                file << "# Wedge factor: " << wedgeFactor << "\n";
                file << "# Fuel density: " << rhoFuel << " kg/m3\n";
                file << "# Air density: " << rhoAir << " kg/m3\n";
                file << "# Inlet area (full): " << inletAreaFull << " m2\n";
                file << "#\n";
                file << "# Time(s)  MassFlow(kg/s)  MomentumFlux_x(N)  PressureForce_x(N)  TotalForce_x(N)\n";
            }
            
            file << mesh.time().value() << "  "
                 << massFlowFull << "  "
                 << momentumFluxFull.x() << "  "
                 << pressureForceFull.x() << "  "
                 << totalForce.x() << "\n";
            
            file.close();
//}}} end code

    return true;
}


bool
Foam::
calcMomentumFluxFunctionObject::write()
{
    if (false)
    {
        printMessage("write calcMomentumFlux");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
calcMomentumFluxFunctionObject::end()
{
    if (false)
    {
        printMessage("end calcMomentumFlux");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

