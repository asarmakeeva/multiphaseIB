/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/
#include "error.H"

#include "IBpressureForceTest.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IBpressureForceTest, 0);

addToRunTimeSelectionTable
(
    forceModel,
    IBpressureForceTest,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
IBpressureForceTest::IBpressureForceTest(
    const dictionary &dict,
    cfdemCloud &sm)
    : forceModel(dict, sm),
      propsDict_(dict.subDict(typeName + "Props")),
      voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),             //mod by alice
      voidfractions_(sm.mesh().lookupObject<volScalarField>(voidfractionFieldName_)), //mod by alice
      velFieldName_(propsDict_.lookup("velFieldName")),
      U_(sm.mesh().lookupObject<volVectorField>(velFieldName_)),
      densityFieldName_(dict_.lookupOrDefault<word>("densityFieldName", "rho")),
      rho_(sm.mesh().lookupObject<volScalarField>(densityFieldName_)),
      pressureFieldName_(propsDict_.lookup("pressureFieldName")),
      p_(sm.mesh().lookupObject<volScalarField>(pressureFieldName_))
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("IBpressureForceForce");  //first entry must the be the force
    particleCloud_.probeM().writeHeader();

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    forceSubM(0).setSwitches(SW_TREAT_FORCE_DEM,true); // treatDEM = true
    Info << "accounting for Archimedes only on DEM side!" << endl;

    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

IBpressureForceTest::~IBpressureForceTest() 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void IBpressureForceTest::setForce() const
{
    vector force;

//to treat gravity induced pressue
   volVectorField gradP = fvc::grad(p_);

 //  volVectorField IBDragPerV = rho_*particleCloud_.turbulence().nu() * fvc::laplacian(U_) - fvc::grad(p_);

#include "setupProbeModel.H"

                                for (int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
 
            force=vector::zero;
            for(int subCell=0;subCell<particleCloud_.voidFractionM().cellsPerParticle()[index][0];subCell++)
            {
                label cellI = particleCloud_.cellIDs()[index][subCell];
                if (cellI > -1) // particle found
                {
                   force += -gradP[cellI]*particleCloud_.mesh().V()[cellI]*(1-voidfractions_[cellI]);
	    
                   //force -= H_[cellI] * H_.mesh().V()[cellI];
                   //force += IBDragPerV[cellI] * IBDragPerV.mesh().V()[cellI];
                }
            }

            //Set value fields and write the probe
            if(probeIt_)
            {
                #include "setupProbeModelfields.H"
                vValues.append(force);           //first entry must the be the force
                particleCloud_.probeM().writeProbe(index, sValues, vValues);
            }

            // set force on particle
            // write particle based data to global array
            forceSubM(0).partToArray(index,force,vector::zero);
        //}
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
