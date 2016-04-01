#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "lpmConfig.h"
#include "OutputMessage.h"
#include "Logger.h"
#include "GlobalConstants.h"
#include "LpmXyzVector.hpp"
#include "LpmEuclideanCoords.hpp"
#include "LpmSphereCoords.hpp"
#include "LpmParticles.hpp"
#include "LpmScalarField.hpp"
#include "LpmField2d.hpp"
#include "LpmField3d.hpp"

typedef double ST;
using LpmXyzVector::XyzVector;

int main ( int argc, const char* argv[] ) {
	Logger* log = Logger::Instance(OutputMessage::debugPriority);

	std::stringstream ss;
    std::string s;
    ss << "-- LPM Version " << LPM_VERSION_MAJOR << "." << LPM_VERSION_MINOR << " -- \n ";
    ss << "   Unit Test " << argv[0] << ": covers LpmParticles.hpp, LpmScalarField.hpp, " 
       << "LpmField2d.hpp and LpmField3d.hpp.";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(introMsg);
    
    /*
    TEST 1 : Planar particles, Cartesian geometry
    */
    ss.str(s);
    ss << "TEST 1: Particles and Fields in the plane" ;
    OutputMessage test1Start( ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(test1Start);
    
    const int nn = 101;
    const ST dx = 0.1;
    
    LpmParticles<ST> planeParticles( LpmParticles<ST>::PlanarCartesian, nn*nn );
    planeParticles.PrintStats( "main, after construction.");
    for ( int i = 0; i < nn; ++i ) {
    	const ST xi = -5.0 + i * dx;
    	for ( int j = 0; j < nn/2; ++j ) {
    		const ST yj = -5.0 + j * dx;
    		planeParticles.insert( xi, yj );
    	}
    	for ( int j = nn/2; j < nn; ++j ) {
			const ST yj = -5.0 + j * dx;
			const XyzVector<ST> physX(xi,yj);
			const XyzVector<ST> lagX(physX);
			planeParticles.insert( physX, lagX );    		
    	}
    }   
    planeParticles.PrintStats( "main, after particle insertion.");
    
    std::string matlabFile = "planarParticlesUnitTest.m";
    std::ofstream file( matlabFile.c_str() );
    if ( !file ) {
    	OutputMessage fileErrMsg("ERROR: cannot open planar output file.", OutputMessage::errorPriority,
    		 "lpmParticlesAndFieldTest.cpp::main");
    	log->logMessage(fileErrMsg);
    	return 1;
    }
    
    planeParticles.writePhysCoordsToMatlab( file );
    planeParticles.writeLagCoordsToMatlab( file );
    planeParticles.writeAreaToMatlab( file );
    
    LpmScalarField<ST> gaussHill( planeParticles, "gaussHill", "n_a" );
    for ( int i = 0; i < planeParticles.size(); ++i ) {
    	const XyzVector<ST> physX = planeParticles.physCoordVec(i);
    	gaussHill.insert( std::exp( -physX.x * physX.x - physX.y * physX.y ) );
    }
    gaussHill.PrintStats( "main, after field value insertion");
    gaussHill.writeFieldToMatlab( file );
    
    LpmField2d<ST> rotatingField( nn * nn, "rotation", "n_a" );
    for ( int i = 0; i < planeParticles.size() / 2; ++i ) {
    	const XyzVector<ST> physX = planeParticles.physCoordVec(i);
    	rotatingField.insert( - physX.y, physX.x );
    }
    for ( int i = planeParticles.size()/2; i < planeParticles.size(); ++i ) {
    	const XyzVector<ST> physX = planeParticles.physCoordVec(i);
    	const XyzVector<ST> fieldVec( - physX.y, physX.x );
    	rotatingField.insert( fieldVec );
    }
    rotatingField.PrintStats("main, after field2d value insertion.");
    rotatingField.writeFieldToMatlab( file );
    
    OutputMessage endTestMsg("TEST COMPLETE", OutputMessage::remarkPriority, "main");
    log->logMessage(endTestMsg);
    
    /*
    TEST 2 : 3d particles, Euclidean geometry
    */
	OutputMessage test2start("TEST 2: 3D Particles and Fields", OutputMessage::remarkPriority, "main");
	log->logMessage(test2start);
    
return 0;
};
