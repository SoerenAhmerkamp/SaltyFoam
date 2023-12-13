	/*---------------------------------------------------------------------------*\
2  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    buoyantBoussinesqPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent flow"
        " of incompressible fluids.\n"
        "Uses the Boussinesq approximation."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"
	#include <math.h>
    turbulence->validate();
	//dimensionedScalar Prs("Prs", dimless, laminarTransport);
	//dimensionedScalar Prst("Prst", dimless, laminarTransport);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"
            #include "SEqn.H"
            #include "NutEqn.H"

	    //rhok = 1.0 - beta*(T - TRef);
		const scalarField& Tvalues =  T.internalField();
		const scalarField& Svalues = S.internalField();
		
		
			forAll(rhok, i)
			{
				rhok[i] = pow(10,3) * 
		(0.5 * (4.032219*0.5 + 0.115313*((2*Svalues[i]-150)/150) + 3.26*pow(10,-4) * (2*pow((2*Svalues[i]-150)/150,2) - 1))
		+ (2*Tvalues[i]-200)/160 * (-0.108199*0.5 + 1.571*pow(10,-3)*((2*Svalues[i]-150)/150) - 4.23*pow(10,-4) * (2*pow((2*Svalues[i]-150)/150,2) - 1))
		+ (2*pow((2*Tvalues[i]-200)/160,2) - 1) * (-0.012247*0.5 + 1.74*pow(10,-3)*((2*Svalues[i]-150)/150) - 9*pow(10,-6) * (2*pow((2*Svalues[i]-150)/150,2) - 1))
		+ (4*pow((2*Tvalues[i]-200)/160,3) - 3*(2*Tvalues[i]-200)/160) * (6.92*pow(10,-4)*0.5 - 8.7*pow(10,-5)*((2*Svalues[i]-150)/150) - 5.3*pow(10,-5) * (2*pow((2*Svalues[i]-150)/150,2) - 1)));
			}
	
		
		// --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
