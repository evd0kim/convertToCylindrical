/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
    convertToCylindrical

Description
    Converts a velocity field from Cartesian coordinates to cylindrical coordinates

Author
     Bryan Lewis, Penn State University
     Assembled from forum threads by Hrvoje Jasak and Hakkan Nilsson
     Ilya Evdokimov, Institute for System Programming of the Russian Academy of Sciences

Usage
     After the simulation has completed, run this application to convert the velocity field to
     cylindrical coordinates (r,theta,z)

     Velocity field must be titled "U"

     The utility uses dynamicMeshDict to read rotation axis and center point automatically

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    argList::addBoolOption
    (
        "unitVectors",
        "save unit vectors of the cylindrical CS"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    //looking for additional utility parameter for save option
    const bool saveUnitVectors = args.optionFound("unitVectors");

	  Info<< "Reading dynamic mesh properties\n" << endl;

    IOdictionary rotationProperties
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    );

	vector axisVec
	(
	    rotationProperties.subDict("solidBodyMotionFvMeshCoeffs").subDict("rotatingMotionCoeffs").lookup("axis")
	);

  point rotCenter
  (
      rotationProperties.subDict("solidBodyMotionFvMeshCoeffs").subDict("rotatingMotionCoeffs").lookup("origin")
  );

	//vector Omega = omega*axisVec;
  const vector globalX(1, 0, 0);
  const vector globalY(0, 1, 0);
  const vector globalZ(0, 0, 1);

  vector dirVec = axisVec;

  if (axisVec == globalX){ dirVec = globalY;};
  if (axisVec == globalY){ dirVec = globalZ;};
  if (axisVec == globalZ){ dirVec = globalX;};

  forAll(timeDirs, timeI)
  {
        runTime.setTime(timeDirs[timeI], timeI);

        Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;

        mesh.readUpdate();
        // defining cylindrical coordinate system

        Info<< "    Creating cylindrical system (r, teta, z)" << endl;

        cylindricalCS cyl
        (
        		"cylindricalCS",
        		rotCenter,      //center
        		axisVec,        //axis
        		dirVec,         //direction
        		false
        ); //false does not play (degree-radians switch)

        //Create unit vectors in the centers of the volumes
        //radial direction unit vector
        volVectorField cRad
        (
              IOobject
              (
                  "cRad",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              (mesh.C() - cmptMultiply(mesh.C(), axisVec))/mag(mesh.C() - cmptMultiply(mesh.C(), axisVec))
  	    );
        //tangential direction unit vector
        volVectorField cTheta
        (
              IOobject
              (
                  "cTheta",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              (axisVec^cRad)/mag(axisVec^cRad)
        );
        //checking save flag
        if (saveUnitVectors)
        {
            Info<< "    Saving unit vectors cRad and cTheta" << endl;
            cRad.write();
            cTheta.write();
        }

        //Set up U
        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (Uheader.headerOk())   // Check that U exists
        {
            mesh.readUpdate();

            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

      	    // Set up Ucyl
      	    volVectorField Ucyl
            (
                  IOobject
                  (
      		           "Ucyl",
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ
                     //IOobject::AUTO_WRITE
                  ),
      		        mesh,
      		        dimensionedVector
      		        (
      		           "Ucyl",
      		           dimensionSet(0,1,-1,0,0,0,0),
      		           vector::zero
      		        )
      	   );
      	    // transformation of velocity field U from cartesian -> cylindrical
      	    Info<< "    Converting U\n" << endl;

            //Transform cc to cylindrical CS
            Ucyl.replace(vector::X, (U&(cRad) ));   //Ur
            Ucyl.replace(vector::Y, (U&(cTheta) ));	//Uteta
            Ucyl.replace(vector::Z, (U&axisVec)); 	//Uz

            Ucyl.write();

	      }
      	else
      	{
        	Info<< "\n    Failed! No exisiting U field\n" << endl;
      	}
    }
    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
