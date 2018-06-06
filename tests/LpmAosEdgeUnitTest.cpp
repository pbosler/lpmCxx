#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAosParticle.hpp"
#include "LpmAosParticleSet.hpp"
#include "LpmAosEdge.hpp"
#include "LpmAosTypes.hpp"
#include "LpmAosEdgeFactory.hpp"
#include "LpmAosEdgeSet.hpp"


using namespace Lpm;

int main (int argc, char* argv[]) {

    typedef std::array<index_type,2> arr_t;

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "Lpm Edge Unit Test" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify Edge and its subclasses' basic functionality." << std::endl;
        ss << "\t 2. Verify EdgeSet basic functionality." << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }
    
    // Edge
    std::shared_ptr<ParticleFactory<3>> sphereFactory(new SWEParticleFactory<3>());
    std::shared_ptr<EdgeFactory<3>> linearEdgeFactory(new LinearEdgeFactory<3>());
    std::shared_ptr<EdgeFactory<3>> quadraticEdgeFactory(new QuadraticEdgeFactory<3>());
    std::shared_ptr<EdgeFactory<3>> cubicEdgeFactory(new CubicEdgeFactory<3>());
    
    ParticleSet<3> cubedSphere(sphereFactory);
    cubedSphere.initFromParticleSetFile("cubed_sphere/cubedSphere_1.vtk");
    std::cout << cubedSphere.infoString() << std::endl;

    std::unique_ptr<Edge<3>> linearEdge = linearEdgeFactory->createEdge(0, 1, 0, 1);
    std::cout << linearEdge->infoString() << std::endl;
    
    std::unique_ptr<Edge<3>> quadEdge = quadraticEdgeFactory->createEdge(0,1,2,3, arr_t({4,-1}));
    std::cout << "quadEdge info: " << std::endl;
    std::cout << quadEdge->infoString() << std::endl;
    
    std::unique_ptr<Edge<3>> cubicEdge = cubicEdgeFactory->createEdge(1,2,0,2, arr_t({1,2}));
    std::cout << "cubic edge info: " << std::endl;
    std::cout << cubicEdge->infoString() << std::endl;

    Vec<3> midpt = linearEdge->midpoint(cubedSphere);
    std::cout << "\torig = " << cubedSphere.physCrd(linearEdge->orig()) << std::endl;
    std::cout << "\tdest = " << linearEdge->destCrd(cubedSphere) << std::endl;
    std::cout << "\tonBoundary = " << (linearEdge->onBoundary() ? "true" : "false") << std::endl;
    std::cout << "\tmidpt= " << midpt << std::endl;
    std::cout << "\tedge_vector = " << linearEdge->edgeVector(cubedSphere) << std::endl;
    std::cout << "\tsphLength = " << linearEdge->sphLength(cubedSphere) << std::endl;
    std::cout << "\teuclength = " << linearEdge->eucLength(cubedSphere) << std::endl;

    // EdgeSet
    EdgeSet<3> edges(linearEdgeFactory, SPHERICAL_SURFACE_GEOMETRY);
    
    std::cout << edges.infoString() << std::endl;
    
return 0;
}
