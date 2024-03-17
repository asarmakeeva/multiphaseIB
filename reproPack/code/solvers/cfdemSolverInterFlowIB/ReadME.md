just reminder what I added to InterFlow solver, 
IsoFoam should be installed, and could be donwloaded from github https://github.com/isoAdvector/isoAdvector

to add paricles into openFOAM solver
`
#include "cfdemCloudIB.H"
#include "implicitCouple.H"

#include "averagingModel.H"
#include "regionModel.H"
#include "voidFractionModel.H"

inside main loop 

cfdemCloudIB particleCloud(mesh);
 then after     while (runTime.loop()) add

        //=== dyM ===================
        interFace = mag(mesh.lookupObject<volScalarField>("voidfractionNext"));
        interFace -= mag(alpha1 * alpha2);
        mesh.update(); //dyM

        // do particle stuff

        Info << "- evolve()" << endl;
        particleCloud.evolve();
        volScalarField voidfractionNext=mesh.lookupObject<volScalarField>("voidfractionNext");`

and inside pimple loop at the end

`
particleCloud.calcCorrectionTerm(U,voidfractionNext,H);
Info << "particleCloud.calcVelocityCorrection() " << endl;
`
also you need to add additional fields in createFields.H

`
volVectorField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
     mesh
);


dimensionedScalar partDens("0", dimensionSet(1, -3, 0, 0, 0), 750);

    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field phiIB\n" << endl;
    volScalarField phiIB
    (
        IOobject
        (
            "phiIB",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field phiIB\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field interFace\n" << endl;
    volScalarField interFace
    (
        IOobject
        (
            "interFace",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    );
`
in alphaEqnSubCycle.H at end add

`
rho == alpha1*rho1 + alpha2*rho2;
rho += (1.0-voidfractionNext)*partDens;
`
add in UEqn + H



in option include

-I$(CFDEM_OFVERSION_DIR) \
-I$(CFDEM_OFVERSION_DIR) \

and include links
-L$(CFDEM_LIB_DIR)\
-l$(CFDEM_LIB_NAME) \
$(CFDEM_ADD_LIB_PATHS) \
$(CFDEM_ADD_LIBS)


after in the folder where main solver code run wmake command






