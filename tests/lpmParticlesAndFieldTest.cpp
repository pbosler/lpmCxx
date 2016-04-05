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

ST sphHarm54( const ST lat, const ST lon) {
	static const ST amp = (3.0 / 16.0 ) * std::sqrt(385.0 / (2.0 * GlobalConstants::Instance()->Pi()));
	return amp * std::cos( 4.0 * lon ) * std::pow( std::sin(lat), 4) * std::cos(lat);
}

int main ( int argc, const char* argv[] ) {
	Logger* log = Logger::Instance(OutputMessage::debugPriority);

	std::stringstream ss;
    std::string s;
    ss << "-- LPM Version " << LPM_VERSION_MAJOR << "." << LPM_VERSION_MINOR << " -- \n ";
    ss << "   Unit Test " << argv[0] << ": covers LpmParticles.hpp, LpmScalarField.hpp, " 
       << "LpmField2d.hpp and LpmField3d.hpp.";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(introMsg);
    
    std::string matlabFile = "particlesAndFieldsUnitTest.m";
	std::ofstream file( matlabFile.c_str() );
    /*
    TEST 1 : Planar particles, Cartesian geometry
    */
    OutputMessage endTestMsg("TEST COMPLETE", OutputMessage::remarkPriority, "main");
    {
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
	
		ss.str(s);
		ss << "distance between particle 0 and particle " << planeParticles.size() - 1 << " should be " 
		   << 10.0 * std::sqrt(2.0) << "; computed distance is " << planeParticles.distance(0, nn * nn - 1 );
		OutputMessage statusMsg(ss.str(), OutputMessage::remarkPriority, "main, plane particles Test 1");
		log->logMessage(statusMsg);
	
		
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
	
		log->logMessage(endTestMsg);
    
    }
    /*
    TEST 2 : 3d particles, Euclidean geometry
    */
    {
		OutputMessage test2start("TEST 2: 3D Particles and Fields", OutputMessage::remarkPriority, "main");
		log->logMessage(test2start);
	
		log->logMessage(endTestMsg);
	}
	/* 
	TEST 3 : Spherical surface
	*/
	{
		OutputMessage test3start("TEST 3: Spherical surface particles and fields", OutputMessage::remarkPriority,
			 "main");
		log->logMessage(test3start);
	
		const ST PI = GlobalConstants::Instance()->Pi();
		const ST deg2rad = GlobalConstants::Instance()->Deg2Rad();
		const int nLat = 91;
		const int nLon = 180;
		const ST dLam = 360.0/nLon * deg2rad;
	
		LpmParticles<ST> sphereParticles( LpmParticles<ST>::SphericalSurface, nLon * nLat );
		for ( int j = 0; j < nLon; ++j ) {
			const ST lon = j * dLam;
			for ( int i = 0; i < nLat; ++i ) {
				const ST lat = -0.5 * PI + i * dLam;
				const XyzVector<ST> sphX( std::cos(lon) * std::cos(lat), std::sin(lon) * std::cos(lat), std::sin(lat));
			
				sphereParticles.insert( sphX, sphX );
			}
		}
		sphereParticles.PrintStats("main, after sphere insertion");
		sphereParticles.writePhysCoordsToMatlab( file );
	
		ss.str(s);
		ss << "distance between particle 0 and particle " << sphereParticles.size() - 1 << " should be " << PI 
		   << "; computed distance is " << sphereParticles.distance(0, nLon * nLat - 1);
		OutputMessage statusMsg( ss.str(), OutputMessage::remarkPriority, "main, sphere particles test.");
		log->logMessage(statusMsg);
	
		LpmScalarField<ST> rhWave( sphereParticles, "rhWave", "n_a");
		for ( int i = 0; i < sphereParticles.size(); ++i ) {
			rhWave.insert( sphHarm54( sphereParticles.Latitude(i), sphereParticles.Longitude(i) ) );
		}
		rhWave.writeFieldToMatlab(file);
		log->logMessage(endTestMsg);
    }
return 0;
};
