#include "LpmCoords.h"
#include "LpmEuclideanCoords.h"
#include "LpmXyzVector.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmSphericalCoords.h"
#include <memory>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include <iostream>
#include <sstream>

using namespace Lpm;

int main (int argc, char* argv[] ) {
    std::unique_ptr<Logger> log = std::unique_ptr<Logger>(new Logger(OutputMessage::debugPriority, "Coords_unitTest_log"));
    
    std::stringstream ss;
    std::string s;
    ss << "-- LPM Version " << LPM_VERSION_MAJOR << "." << LPM_VERSION_MINOR << " -- \n ";
    ss << "   Unit Test " << argv[0] << ": covers LpmCoords.h, LpmEuclideanCoords.h, and LpmSphericalCoords.h.";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(introMsg);
    ss.str(s);
    
    //
    // Test Euclidean coordinates
    //
    ss.str(s);
    ss << "TEST 1 : EUCLIDEAN COORDS";
    OutputMessage test1Start( ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(test1Start);
    
    EuclideanCoords ec(13);
    
    ec.insert( -1.0, 1.0 ); //0
    ec.insert( -1.0, 0.0 ); //1
    ec.insert( -1.0,-1.0 ); //2
	ec.insert( 0.0, -1.0 ); //3
	ec.insert( -100.0, 100.0 ); //4: 1.0		-1.0 
	ec.insert(1.0, 0.0 ); // 5
	ec.insert(1.0, 1.0 ); //6
	ec.insert(0.0, 1.0 ); //7
	ec.insert(0.0, 0.0 ); //8
	ec.insert(-0.5, 0.5 ); //9
	
	ec.replace(4, 1.0, -1.0 );
	
	XyzVector vec1(-0.5, -0.5 );
	XyzVector vec2(0.5,	-0.5 );
	XyzVector vec3(0.5,  0.5 );
	
	ec.insert(vec1); //10
	ec.insert(vec2); //11
	ec.insert(vec3); //12
	
	std::string crdStr = ec.listAllCoords();
	std::cout << crdStr << std::endl;
    
    std::cout << "dist between coord 0 and coord 4 should be " << 2.0 * std::sqrt(2.0) 
			  << "; computed distance is " << ec.distance(0,4) << std::endl;
	std::cout << "midpoint between coord 2 and coord 6 should be the origin; computed midpoint is " 
			  << ec.midpoint(2,6) << std::endl;
	std::vector<index_type> inds;
	inds.push_back(8);
	inds.push_back(1);
	inds.push_back(3);
	inds.push_back(2);
	std::cout << "centroid should be (-0.5, -0.5); computed centroid is " << ec.centroid(inds) << std::endl;	  
	std::cout << "area of triangle 2, 6, 0 should be 2.0; computed area is " << ec.triArea(2,6,0) << std::endl;
	
	OutputMessage endTest( "END TEST", OutputMessage::remarkPriority, "main");
    log->logMessage(endTest);
    
    
    OutputMessage test2start( "TEST 2 : SPHERE COORDS ", OutputMessage::remarkPriority, "main");
    log->logMessage(test2start);
    
    SphericalCoords sc(14);
    sc.insert( 0.577350269189626,	-0.577350269189626,	0.577350269189626);//0
	sc.insert( 0.577350269189626,	-0.577350269189626,	-0.577350269189626);//1
	sc.insert( 0.577350269189626,	0.577350269189626,	-0.577350269189626);//2
	sc.insert( 0.577350269189626,	0.577350269189626,	0.577350269189626);//3
	sc.insert( -0.577350269189626,	0.577350269189626,	-0.577350269189626);//4
	sc.insert( -0.577350269189626,	0.577350269189626,	0.577350269189626);//5
	sc.insert( -0.577350269189626,	-0.577350269189626,	-0.577350269189626);//6
	sc.insert( -0.577350269189626,	-0.577350269189626,	0.577350269189626);//7
	sc.insert( 1.0, 0.0, 0.0);//8
	sc.insert( 0.0, 1.0, 0.0);//9
	sc.insert( -1.0, 0.0, 0.0);//10
	sc.insert( 0.0, -1.0, 0.0);//11
	sc.insert( 0.0, 0.0, 1.0);//12
	sc.insert( 0.0, 0.0, -1.0);//13
	
	
	for ( int i = 0; i < sc.n(); ++i )
		std::cout << "coord " << i << " has magnitude = " << sc.magnitude(i) << std::endl;
	
	inds[0] = 0;
	inds[1] = 3;
	inds[2] = 5;
	inds[3] = 7;
	std::cout << "centroid of indices 0, 3, 5, 7 should be (0,0,1); computed centroid is " 
			  << sc.centroid(inds) << std::endl;
	
	scalar_type surfArea = 0.0;
	// face 1 area
	surfArea += sc.triArea(0,8,1);
	surfArea += sc.triArea(1,8,2);
	surfArea += sc.triArea(2,8,3);
	surfArea += sc.triArea(3,8,0);
	// face 2
	surfArea += sc.triArea(3,9,2);
	surfArea += sc.triArea(2,9,4);
	surfArea += sc.triArea(4,9,5);
	surfArea += sc.triArea(5,9,3);
	// face 3
	surfArea += sc.triArea(5,10,4);
	surfArea += sc.triArea(4,10,6);
	surfArea += sc.triArea(6,10,7);
	surfArea += sc.triArea(7,10,5);
	// face 4
	surfArea += sc.triArea(7,11,6);
	surfArea += sc.triArea(6,11,1);
	surfArea += sc.triArea(1,11,0);
	surfArea += sc.triArea(0,11,7);
	// face 5
	surfArea += sc.triArea(7,12,0);
	surfArea += sc.triArea(0,12,3);
	surfArea += sc.triArea(3,12,5);
	surfArea += sc.triArea(5,12,7);
	// face 6
	surfArea += sc.triArea(1,13,6);
	surfArea += sc.triArea(6,13,4);
	surfArea += sc.triArea(4,13,2);
	surfArea += sc.triArea(2,13,1);
	
 	std::cout << "sphere surface area should be " << 4.0 * PI << "; computed area is "
 			  << surfArea << std::endl;
 			  
 	std::cout << "midpoint of coord 0 and coord 2 should be (1,0,0); computed midpoint is " 
 			  << sc.midpoint(0,2) << std::endl;		  
	
	log->logMessage(endTest);
return 0;
}