#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include "lpmConfig.h"
#include "OutputMessage.h"
#include "Logger.h"
#include "LpmXyzVector.hpp"
#include "LpmEuclideanCoords.hpp"
#include "LpmSphereCoords.hpp"
#include "GlobalConstants.h"

typedef double ST;
using LpmXyzVector::XyzVector;

int main ( int argc, const char* argv[] ) {
	Logger* log = Logger::Instance(OutputMessage::debugPriority);

	std::stringstream ss;
    std::string s;
    ss << "-- LPM Version " << LPM_VERSION_MAJOR << "." << LPM_VERSION_MINOR << " -- \n ";
    ss << "   Unit Test " << argv[0] << ": covers LpmCoords.hpp and LpmSphereCoords.hpp.";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(introMsg);
    
    //
    // Test standard coordinates in the plane
    //
    ss.str(s);
    ss << "TEST 1 : PLANE COORDS";
    OutputMessage test1Start( ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(test1Start);
    
    LpmEuclideanCoords<ST> planeCoords( 2, 13 );
    
    planeCoords.insert( -1.0, 1.0 ); //0
    planeCoords.insert( -1.0, 0.0 ); //1
    planeCoords.insert( -1.0,-1.0 ); //2
	planeCoords.insert( 0.0, -1.0 ); //3
	planeCoords.insert( -100.0, 100.0 ); //4: 1.0		-1.0 
	planeCoords.insert(1.0, 0.0 ); // 5
	planeCoords.insert(1.0, 1.0 ); //6
	planeCoords.insert(0.0, 1.0 ); //7
	planeCoords.insert(0.0, 0.0 ); //8
	planeCoords.insert(-0.5, 0.5 ); //9
	
	XyzVector<ST> vec1(-0.5, -0.5 );
	XyzVector<ST> vec2(0.5,	-0.5 );
	XyzVector<ST> vec3(0.5,  0.5 );
	
	planeCoords.insert(vec1); //10
	planeCoords.insert(vec2); //11
	planeCoords.insert(vec3); //12
	
	planeCoords.replace(4, 1.0, -1.0 );
	
	std::cout << planeCoords << std::endl;
	
	std::cout << "dist between coord 0 and coord 4 should be " << 2.0 * std::sqrt(2.0) 
			  << "; computed distance is " << planeCoords.distance(0,4) << std::endl;
	std::cout << "midpoint between coord 2 and coord 6 should be the origin; computed midpoint is " 
			  << planeCoords.midpoint(2,6) << std::endl;
	std::vector<size_t> inds;
	inds.push_back(8);
	inds.push_back(1);
	inds.push_back(3);
	inds.push_back(2);
	std::cout << "centroid should be (-0.5, -0.5); computed centroid is " << planeCoords.centroid(inds) << std::endl;	  
	std::cout << "area of triangle 2, 6, 0 should be 2.0; computed area is " << planeCoords.triArea(2,6,0) << std::endl;
    
    OutputMessage endTest( "END TEST", OutputMessage::remarkPriority, "main");
    log->logMessage(endTest);
    
    OutputMessage test2start( "TEST 2 : SPHERE COORDS ", OutputMessage::remarkPriority, "main");
    log->logMessage(test2start);
    
    //
    //	Test spherical coordinates
    //
    LpmSphereCoords<ST> sphereCoords( 14 );
	sphereCoords.insert( 0.577350269189626,	-0.577350269189626,	0.577350269189626);//0
	sphereCoords.insert( 0.577350269189626,	-0.577350269189626,	-0.577350269189626);//1
	sphereCoords.insert( 0.577350269189626,	0.577350269189626,	-0.577350269189626);//2
	sphereCoords.insert( 0.577350269189626,	0.577350269189626,	0.577350269189626);//3
	sphereCoords.insert( -0.577350269189626,	0.577350269189626,	-0.577350269189626);//4
	sphereCoords.insert( -0.577350269189626,	0.577350269189626,	0.577350269189626);//5
	sphereCoords.insert( -0.577350269189626,	-0.577350269189626,	-0.577350269189626);//6
	sphereCoords.insert( -0.577350269189626,	-0.577350269189626,	0.577350269189626);//7
	sphereCoords.insert( 1.0, 0.0, 0.0);//8
	sphereCoords.insert( 0.0, 1.0, 0.0);//9
	sphereCoords.insert( -1.0, 0.0, 0.0);//10
	sphereCoords.insert( 0.0, -1.0, 0.0);//11
	sphereCoords.insert( 0.0, 0.0, 1.0);//12
	sphereCoords.insert( 0.0, 0.0, -1.0);//13
	
	std::cout << sphereCoords << std::endl;
	for ( int i = 0; i < sphereCoords.size(); ++i )
		std::cout << "coord " << i << " has magnitude = " << sphereCoords.magnitude(i) << std::endl;
	
	inds[0] = 0;
	inds[1] = 3;
	inds[2] = 5;
	inds[3] = 7;
	std::cout << "centroid of indices 0, 3, 5, 7 should be (0,0,1); computed centroid is " 
			  << sphereCoords.centroid(inds) << std::endl;
			  
	static const ST PI = GlobalConstants::Instance()->Pi();
	ST surfArea = 0.0;
	// face 1 area
	surfArea += sphereCoords.triArea(0,8,1);
	surfArea += sphereCoords.triArea(1,8,2);
	surfArea += sphereCoords.triArea(2,8,3);
	surfArea += sphereCoords.triArea(3,8,0);
	// face 2
	surfArea += sphereCoords.triArea(3,9,2);
	surfArea += sphereCoords.triArea(2,9,4);
	surfArea += sphereCoords.triArea(4,9,5);
	surfArea += sphereCoords.triArea(5,9,3);
	// face 3
	surfArea += sphereCoords.triArea(5,10,4);
	surfArea += sphereCoords.triArea(4,10,6);
	surfArea += sphereCoords.triArea(6,10,7);
	surfArea += sphereCoords.triArea(7,10,5);
	// face 4
	surfArea += sphereCoords.triArea(7,11,6);
	surfArea += sphereCoords.triArea(6,11,1);
	surfArea += sphereCoords.triArea(1,11,0);
	surfArea += sphereCoords.triArea(0,11,7);
	// face 5
	surfArea += sphereCoords.triArea(7,12,0);
	surfArea += sphereCoords.triArea(0,12,3);
	surfArea += sphereCoords.triArea(3,12,5);
	surfArea += sphereCoords.triArea(5,12,7);
	// face 6
	surfArea += sphereCoords.triArea(1,13,6);
	surfArea += sphereCoords.triArea(6,13,4);
	surfArea += sphereCoords.triArea(4,13,2);
	surfArea += sphereCoords.triArea(2,13,1);
	
 	std::cout << "sphere surface area should be " << 4.0 * PI << "; computed area is "
 			  << surfArea << std::endl;
 			  
 	std::cout << "midpoint of coord 0 and coord 2 should be (1,0,0); computed midpoint is " 
 			  << sphereCoords.midpoint(0,2) << std::endl;		  
	
	log->logMessage(endTest);
return 0;
};
