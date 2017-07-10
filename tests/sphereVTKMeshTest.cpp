#include "GlobalConstants.h"
#include "PolyMesh2d.h"
#include "OutputMessage.h"
#include "Logger.h"
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
	
	for (int i = 1; i < 9; ++i) 
	{
		//
		// construct 2 sphere mesh objects
		//
// 		const int i = 1;
		const int initNest = i;
		const double sphereRadius = 1.0;
	// 	constants->SetEarthRadius( sphereRadius );
	
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
		ss.str("");
// 		ss << "icosTri_" << i << ".vtk";
// 		triMesh.outputToVTK( ss.str() );
 		ss << "icosTri_" << i << ".vtk";
        triMesh.writeActiveParticlesToVTK(ss.str(), "icosahedral triangles");
        std::cout << "write routine finished, back to main.\n";
		ss.str("");
// 		ss << "cubedSphere_" << i << ".vtk";
// 		quadMesh.outputToVTK(ss.str());
        ss << "cubedSphere_" << i << ".vtk";
        quadMesh.writeActiveParticlesToVTK(ss.str(), "spherical quadrilaterals");
        
        std::cout << "write routine finished, back to main.\n";
	}
	
	return 0;
};

