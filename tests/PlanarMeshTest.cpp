/** @file PlanarMeshTester.cpp
	@author Peter Bosler, Sandia National Laboratories
	@brief Unit test for PolyMesh2d objects and implementation.
*/
#include "PolyMesh2d.h"
#include "xyzVector.h"
#include "Field.h"
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
    ss << "		 This test builds planar meshes. \n";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");
    exeLog->logMessage(introMsg);	
	
	const int initNest = 3;
	const double maxR = 1.0;
	
	std::cout << "*** Testing basic hexagonal/tri mesh construction ... " << std::endl;
	
	//
	// Test the basic constructor with a hexagon/triangle mesh
	//
	PolyMesh2d hexMesh(initNest, PolyMesh2d::triHexSeed, maxR, procRank, numProcs );
	LongMessage hexMeshStatus("hex/tri mesh built", OutputMessage::tracePriority, "main", hexMesh.getInfo() );
	exeLog->logMessage(hexMeshStatus);
	
	std::cout << "hexMesh has " << hexMesh.nDividedFaces() << " divided faces." << std::endl;
	

	char filename[56];
	int j = std::sprintf(filename, "hexMesh_nest%d.m", initNest);
	hexMesh.outputToMatlab( filename );
	std::cout << "*** Test complete" << std::endl;

	//
	// test basic constructor with a rectangular/quadrilateral mesh
	//	
	std::cout << "*** Testing basic rectangular/quad mesh construction ... " << std::endl;
	PolyMesh2d quadMesh(initNest, PolyMesh2d::quadRectSeed, maxR, procRank, numProcs);	
	LongMessage quadMeshStatus("rect/quad mesh built", OutputMessage::tracePriority, "main", quadMesh.getInfo() );
	exeLog->logMessage(quadMeshStatus);
	
	std::cout << "quadMesh has " << quadMesh.nDividedFaces() << " divided faces." << std::endl;
	
	j = std::sprintf(filename, "quadMesh_nest%d.m", initNest);
	quadMesh.outputToMatlab( filename);
	
	std::cout << "*** Test complete ... " << std::endl;
	
	//
	// test the point-query algorithm
	//
	std::cout << "*** Testing point query algorithm in both meshes ... " << std::endl;
	const int N = 100;
	double x[N+1];
	double y[N+1];
	for ( int i = 0; i <= N; ++i )
	{
		x[i] = -maxR + 2.0*maxR * i / N;
		y[i] = -maxR + 2.0*maxR * i / N;
	}
	
	int hexlocs[N+1][N+1];
	int quadlocs[N+1][N+1];
	for ( int j =0; j <= N; ++j )
	{
		for ( int i = 0; i <= N; ++i )		
		{
			xyzVector xij( x[i], y[j] );			
		
			hexlocs[i][j] = hexMesh.locatePointInMeshFace( xij );
			quadlocs[i][j] = quadMesh.locatePointInMeshFace( xij );
		}
	}
	std::cout << "... query function returned ... " << std::endl;
	std::string locFile = "locations.m";
	
	std::ofstream hexfile( locFile.c_str() );
	hexfile << "xij = [";

	for ( int i = 0; i < N; ++i )
	{
		hexfile << x[i] << ", ";
	}
	hexfile << x[N] << "];" << std::endl;
	
	hexfile << "hexLoc = [";
	for ( int i = 0; i < N; ++i )
	{
		for ( int j = 0; j < N; ++j )
		{
			hexfile << hexlocs[i][j] << ", ";
		}
		hexfile << hexlocs[i][N] << "; ";
	}
	for ( int j = 0; j < N; ++j)
	{
		hexfile << hexlocs[N][j] << ", ";
	}
	hexfile << hexlocs[N][N] << "]; " << std::endl;

	hexfile << "quadLoc = [";
	for ( int i = 0; i < N; ++i )
	{
		for ( int j = 0; j < N; ++j )
		{
			hexfile << quadlocs[i][j] << ", ";
		}
		hexfile << quadlocs[i][N] << "; ";
	}
	for ( int j = 0; j < N; ++j)
	{
		hexfile << quadlocs[N][j] << ", ";
	}
	hexfile << quadlocs[N][N] << "]; " << std::endl;
	std::cout << "*** Test complete. " << std::endl;

	//
	// test the adjacency functions
	//
	std::cout << "*** Testing mesh adjacency and dual mesh adjacency functions ... " << std::endl;
// 	int hexMeshFaceQuery = hexMesh.nFaces()-1;
// 	int quadMeshFaceQuery = quadMesh.nFaces()-1;
	int hexMeshFaceQuery = 3*hexMesh.nFaces()/4-11;
 	int quadMeshFaceQuery = 3*quadMesh.nFaces()/4-11;
	xyzVector faceXY;
	
	std::string adjacencyFile = "adjacencies.m";
	
	std::ofstream adjFile( adjacencyFile.c_str() );
	
	faceXY = hexMesh.facePosition(hexMeshFaceQuery);
	adjFile << "hexFaceXY = [" << faceXY.x << ", " << faceXY.y << "];"<< std::endl;

	faceXY = quadMesh.facePosition(quadMeshFaceQuery);	
	adjFile << "quadFaceXY = [" << faceXY.x << ", " << faceXY.y << "];" << std::endl;
	
	std::vector<int> hexEdgesAroundFace = hexMesh.ccwEdgesAroundFace( hexMeshFaceQuery );
	std::vector<int> hexAdjFaces = hexMesh.ccwAdjacentFaces( hexMeshFaceQuery );
	std::vector<int> hexVertsAroundFace = hexMesh.ccwVerticesAroundFace( hexMeshFaceQuery );
	int hexVertQuery = hexVertsAroundFace[0];
	std::vector<int> hexFacesAroundVertex = hexMesh.ccwFacesAroundVertex( hexVertQuery );
	
	std::vector<int> quadEdgesAroundFace = quadMesh.ccwEdgesAroundFace( quadMeshFaceQuery );
	std::vector<int> quadAdjFaces = quadMesh.ccwAdjacentFaces( quadMeshFaceQuery );
	std::vector<int> quadVertsAroundFace = quadMesh.ccwVerticesAroundFace( quadMeshFaceQuery );
	int quadVertQuery = quadVertsAroundFace[1];
	std::vector<int> quadFacesAroundVertex = quadMesh.ccwFacesAroundVertex( quadVertQuery );
	
	std::cout << "Querying hexMesh face " << hexMeshFaceQuery << ", and hexMesh vertex " << hexVertQuery << std::endl;
	std::cout << "... face " << hexMeshFaceQuery << " adjFaces = ";
	for ( int i = 0; i < hexAdjFaces.size(); ++i)
		std::cout << hexAdjFaces[i] << ", ";
	std::cout << std::endl;
	
	std::cout << "... faces around vertex " << hexVertQuery << "  = " ;
	for ( int i = 0; i < hexFacesAroundVertex.size(); ++i)
		std::cout << hexFacesAroundVertex[i] << ", ";
	std::cout << std::endl;
	
	xyzVector vertXY = hexMesh.particlePosition(hexVertQuery);
	adjFile << "hexVertXY = [" << vertXY.x << ", " << vertXY.y << "];" << std::endl;
	adjFile << "hexFacesAroundVertXY = [" ;
	for ( int i = 0; i < hexFacesAroundVertex.size()-1; ++i)
	{
		faceXY = hexMesh.facePosition( hexFacesAroundVertex[i] );
		adjFile << faceXY.x << ", " << faceXY.y << "; ";
	}
	faceXY = hexMesh.facePosition( hexFacesAroundVertex[hexFacesAroundVertex.size()-1] );
	adjFile << faceXY.x << ", " << faceXY.y << "];" << std::endl;

	adjFile << "hexAdjFaces = [";
	for ( int i = 0; i < hexAdjFaces.size()-1; ++i )
	{
		faceXY = hexMesh.facePosition( hexAdjFaces[i] );
		adjFile << faceXY.x << ", " << faceXY.y << "; ";
	}
	faceXY = hexMesh.facePosition( hexAdjFaces[hexAdjFaces.size()-1]);
	adjFile << faceXY.x << ", " << faceXY.y << "];" << std::endl;
	
	std::cout << "Querying quadMesh face " << quadMeshFaceQuery << ", and quadMesh vertex " << quadVertQuery << std::endl;
	std::cout << "... face " << quadMeshFaceQuery << " adjFaces = " ;
	for ( int i = 0; i < quadAdjFaces.size(); ++i )
		std::cout << quadAdjFaces[i] << ", ";
	std::cout << std::endl;
	
	vertXY = quadMesh.particlePosition(quadVertQuery );
	adjFile << "quadVertXY = [" << vertXY.x << ", " << vertXY.y << "];" << std::endl;
	adjFile << "quadFacesAroundVertexXY = [";
	for ( int i = 0; i < quadFacesAroundVertex.size() - 1; ++i )
	{
		faceXY = quadMesh.facePosition( quadFacesAroundVertex[i] );
		adjFile << faceXY.x << ", " << faceXY.y << "; " ;
	}
	faceXY = quadMesh.facePosition( quadFacesAroundVertex[quadFacesAroundVertex.size() - 1] );
	adjFile << faceXY.x << ", " << faceXY.y << "]; " << std::endl; 
	
	adjFile << "quadAdjFaces = [";
	for ( int i = 0; i < quadAdjFaces.size()-1; ++i)
	{
		faceXY = quadMesh.facePosition( quadAdjFaces[i] );
		adjFile << faceXY.x << ", " << faceXY.y << "; ";
	}
	faceXY = quadMesh.facePosition( quadAdjFaces[quadAdjFaces.size()-1] );
	adjFile << faceXY.x << ", " << faceXY.y << "]; " << std::endl;
	
	
	std::cout << "*** Test complete. " << std::endl;
return 0;
};