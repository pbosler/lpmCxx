#include "GlobalConstants.h"
#include "PolyMesh2d.h"
#include "OutputMessage.h"
#include "Logger.h"
#include "Field.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "lpmConfig.h"
#include <sstream>

int main ( int argc, char* argv[] )
{
	const int procRank = 0;
	const int numProcs = 1;
	const OutputMessage::priority loggingLevel = OutputMessage::debugPriority;
	
	Logger* exeLog = Logger::Instance(loggingLevel, procRank, numProcs);
	
	std::stringstream ss;
    std::string s;
    ss << "-- LPM Version " << LPM_VERSION_MAJOR << "." << LPM_VERSION_MINOR << " -- \n ";
    ss << "   Unit Test " << argv[0] << " \n";
    ss << "		 This test builds spherical meshes. \n";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");
    exeLog->logMessage(introMsg);	
	
	GlobalConstants* constants = GlobalConstants::Instance();
	
	//
	// construct 2 sphere mesh objects
	//
	const int initNest = 3;
	const double sphereRadius = 1.0;
	constants->SetEarthRadius( sphereRadius );
	
	std::cout << "building icos tri sphere mesh ... " << std::endl;
	PolyMesh2d triMesh( initNest, PolyMesh2d::icosTriSphereSeed, sphereRadius, procRank, numProcs );
	
	std::cout << "icos. tri. mesh built. surf area = " << triMesh.surfaceArea() << std::endl;
	Particles* triParticles = triMesh.getParticles();
	
	//triParticles->printIncidentEdges();
	
	std::cout << "building cubed sphere mesh ... " << std::endl;
	PolyMesh2d quadMesh( initNest, PolyMesh2d::cubedSphereSeed,  sphereRadius, procRank, numProcs );
	std::cout << "cubed sphere built. surf area = " << quadMesh.surfaceArea() << std::endl;
	
	Particles* quadParticles = quadMesh.getParticles();
	//quadParticles->printIncidentEdges();
	
	//
	// output meshes to paraview data files
	//
	triMesh.outputToVTK( "icosTri_3.vtk" );
	quadMesh.outputToVTK("cubedSphere_3.vtk");
	
	
	return 0;
};

