#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <exception>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAosTypes.hpp"

using namespace Lpm::Aos;

using Lpm::index_type;
using Lpm::scalar_type;
using Lpm::Logger;
using Lpm::OutputMessage;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "LPM Array of Structures types unit test" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. verify template classes' basic functions" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }

    Vec<3> v3def;
    std::cout << v3def << std::endl;

    scalar_type xxA[3] = {0.0, 1.0, 2.0};
    scalar_type xxB[3] = {1.0, -1.0, 5.0};

    std::vector<scalar_type> xxC = {-10.0, -100.0, -1000.0};

    Vec<2> v2a(xxA);
    Vec<3> v3a(xxA);
    Vec<3> v3copy(v3a);
    Vec<2> v2fromVec(xxC);

    std::cout << "v2a = " << v2a << std::endl;
    std::cout << "|v2a| = " << v2a.mag() << std::endl;
    std::cout << "v3a    = " << v3a << std::endl;
    std::cout << "v3copy = " << v3copy << std::endl;
    std::cout << "|v3a| = " << v3a.mag() << std::endl;

    std::cout << "v2fromVec = " << v2fromVec << std::endl;

    Vec<2> v2b(xxB);
    scalar_type cp2 = v2a.crossProdComp3(v2b);
    std::cout << "v2b = " << v2b << std::endl;
    std::cout << "v2a cross v2b = " << cp2 << std::endl;
    std::cout << "v2b - v2a = " << v2b - v2a << std::endl;

    Vec<3> v3b(xxB);
    Vec<3> vsum = v3a + v3b;
    Vec<3> cp3 = v3a.crossProd(v3b);
    std::cout << "v3b = " << v3b << std::endl;
    std::cout << "v3a + v3b = " << vsum << std::endl;
    std::cout << "v3a cross v3b = " << v3a.crossProd(v3b) << std::endl;

    std::cout << "v3a.midpoint(v3b) = " << v3a.midpoint(v3b) << std::endl;
    std::cout << "v3a.dotProd(v3b) = " << v3a.dotProd(v3b) << std::endl;
    std::cout << "v3a / v3c = " << v3a / Vec<3>(xxC) << std::endl;
    std::cout << "v3a * v3c = " << v3a * Vec<3>(xxC) << std::endl;
    std::cout << "v3a = v3copy : " << (v3a == v3copy ? "True" : "False") << std::endl;

    Vec<3> v3assign = v3a;

    std::cout << "v3assign = " << v3assign << std::endl;

    std::vector<Vec<2>> vec2s = {v2a, v2b, v2fromVec};
    std::cout << "v2a.midpoint(v2b) = " << v2a.midpoint(v2b) << std::endl;
    std::cout << "plane tri area = " << triArea(vec2s) << std::endl;
    std::cout << "v2b.dist(v2a) = " << v2b.dist(v2a) <<std::endl;
    std::cout << "plane barycenter = " << baryCenter(vec2s) << std::endl;

    Vec<3> v3spha = v3a.normalize();
    Vec<3> v3sphb = v3b.normalize();
    std::cout << "v3a.normalize() " << v3a.normalize() << std::endl;
    std::cout << "v3a.normalize().magSq() = " << v3a.normalize().magSq() << std::endl;
    std::cout << "v3spha.longitude(), latitude() = (" << v3spha.longitude() << ", " << v3spha.latitude() << ")" << std::endl;
    std::cout << "v3spha.sphereMidpoint(v3sphb), magSq = " << v3spha.sphereMidpoint(v3sphb) << ", " << v3spha.sphereMidpoint(v3sphb).magSq() << std::endl;
    std::cout << "v3spha.sphereDist(v3sphB) = " << v3spha.sphereDist(v3sphb) << std::endl;
    std::cout << "v3spha.dist(v3sphB) = " << v3spha.dist(v3sphb) << std::endl;

    std::vector<Vec<3>> vecs = {v3a, v3b, Vec<3>(xxC)};
    std::cout << "barycenter = " << baryCenter(vecs) << std::endl;
    std::cout << "triArea = " << triArea(vecs) << std::endl;

    for (int i=0; i<3; ++i)
        vecs[i].normalizeInPlace();

    std::cout << "sphereBaryCenter = " << sphereBaryCenter(vecs) << std::endl;
    std::cout << "sphereTriArea = " << sphereTriArea(vecs) << std::endl;
return 0;
}

