#include <iostream>
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
return 0;
};
