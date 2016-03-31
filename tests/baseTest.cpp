//
//  main.cpp
//  TestingCode
//
//  Created by Peter Bosler on 11/2/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

/** @file baseTest.cpp
	@author Peter Bosler, Sandia National Laboratories
	@brief Unit test for foundation classes OutputMessage, Logger, and XyzVector<ST>.
*/


#include <iostream>
#include <vector>
#include "OutputMessage.h"
#include "Logger.h"
#include "LpmXyzVector.hpp"
#include "lpmConfig.h"
#include <sstream>

using std::cout;
using std::endl;
using LpmXyzVector::XyzVector;

typedef double ST;

int main(int argc, const char * argv[]) {
    
    Logger* log = Logger::Instance(OutputMessage::debugPriority);
    
    std::stringstream ss;
    std::string s;
    ss << "-- LPM Version " << LPM_VERSION_MAJOR << "." << LPM_VERSION_MINOR << " -- \n ";
    ss << "   Unit Test " << argv[0] << ": covers OutputMessage.h, Logger.h, LpmXyzVector.hpp, and Cmake Config.";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");
    log->logMessage(introMsg);
    
    
    OutputMessage statusMsg("Testing XyzVector<ST> class :", OutputMessage::remarkPriority, "main");
    log->logMessage(statusMsg);
    
    XyzVector<ST> vec3( 0.0, 0.0, 4.0);
    cout << " threeD vector : \n" << vec3 ;
    cout << " has magnitude = " << vec3.magnitude() << endl;
    XyzVector<ST> vec32( 1.0, 0.0, 0.0);
    XyzVector<ST> cp = vec3.crossProduct(vec32);
    cout << " cross product vector : \n" << cp ;
    cout << " has magnitude = " << cp.magnitude() << endl;
    cout << "\n";
      
    cout << "testing copy : " << endl;
    XyzVector<ST> vec3copy(vec3);
    cout << vec3copy << endl;
    
    cout <<  "testing assignment : " << endl;
    vec3 = vec32;
    cout << vec3 << endl;
    
    XyzVector<ST> a(0.0, 0.0);
    XyzVector<ST> b(1.0, 0.0);
    XyzVector<ST> c(1.0, 1.0);
    XyzVector<ST> d(0.0, 1.0);
    
    std::vector< XyzVector<ST> > tri(3);
    tri[0] = a;
    tri[1] = b;
    tri[2] = c;
    std::vector< XyzVector<ST> > quad(4);
    quad[0] = a;
    quad[1] = b;
    quad[2] = c;
    quad[3] = d;
    
    cout << "triangle area should be 0.5; computed area is " << triArea(a, b, c ) << " .\n";
    
    cout << "quad centroid should be (0.5, 0.5); computed centroid is " << centroid(quad) << ".\n";
    cout << "test += operator is part of centroid function.\n";
    
    cout << "testing + operator : " << endl;
    XyzVector<ST> oneone = b + d;
    cout << "sum should be (1,1); sum is " << oneone << ".\n";
    
    cout << "testing - operator : " << endl;
    XyzVector<ST> orig = oneone - b - d;
    cout << "difference should be (0,0); difference is " << orig << ".\n";
    
    cout << "test == operator : " << endl;
    bool ans = (orig == a);
    cout << "answer should be true; answer is " << ans << endl;
    
    statusMsg = OutputMessage("... end test.", OutputMessage::remarkPriority, "main");
    log->logMessage(statusMsg);
    
    return 0;
}
