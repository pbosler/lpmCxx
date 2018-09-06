#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAosParticle.hpp"
#include "LpmAosParticleFactory.hpp"
#include "LpmAosParticleSet.hpp"
#include "LpmAosEdge.hpp"
#include "LpmAosTypes.hpp"
#include "LpmAosEdgeFactory.hpp"
#include "LpmAosEdgeSet.hpp"


using namespace Lpm::Aos;

using Lpm::index_type;
using Lpm::scalar_type;
using Lpm::Logger;
using Lpm::OutputMessage;
using Lpm::GeometryType::SPHERICAL_SURFACE_GEOMETRY;

int main (int argc, char* argv[]) {

    typedef std::vector<index_type> arr_t;

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
    
    ParticleSet<3> cubedSphere(sphereFactory, 100);
    cubedSphere.initFromParticleSetFile("cubed_sphere/cubedSphere_1.vtk");
    std::cout << cubedSphere.infoString() << std::endl;

    std::unique_ptr<Edge<3>> linearEdge = linearEdgeFactory->createEdge(0, 1, 0, 1);
    std::cout << linearEdge->infoString() << std::endl;
    
    KidEdgeArrays<3> linearkids = linearEdge->divide(cubedSphere, 1.0, SPHERICAL_SURFACE_GEOMETRY);
    std::cout << "KidEdgeArrays<3> info:" << std::endl;
    std::cout << linearkids.infoString() << std::endl;
    
    std::unique_ptr<Edge<3>> quadEdge = quadraticEdgeFactory->createEdge(0,1,2,3, arr_t({4,-1}));
    std::cout << "quadEdge info: " << std::endl;
    std::cout << quadEdge->infoString() << std::endl;
    KidEdgeArrays<3> quadkids = quadEdge->divide(cubedSphere, 1.0, SPHERICAL_SURFACE_GEOMETRY);
    std::cout << quadkids.infoString() << std::endl;
    
    std::unique_ptr<Edge<3>> cubicEdge = cubicEdgeFactory->createEdge(1,2,0,2, arr_t({1,2}));
    std::cout << "cubic edge info: " << std::endl;
    std::cout << cubicEdge->infoString() << std::endl;
    
    KidEdgeArrays<3> cubicKids = cubicEdge->divide(cubedSphere, 1.0, SPHERICAL_SURFACE_GEOMETRY);
    std::cout << cubicKids.infoString() << std::endl;

	std::cout << "LinearEdge<3> info:" << std::endl;
    Vec<3> midpt = linearEdge->midpoint(cubedSphere);
    std::cout << "\torig = " << cubedSphere.physCrd(linearEdge->orig()) << std::endl;
    std::cout << "\tdest = " << linearEdge->destCrd(cubedSphere) << std::endl;
    std::cout << "\tonBoundary = " << (linearEdge->onBoundary() ? "true" : "false") << std::endl;
    std::cout << "\tmidpt= " << midpt << std::endl;
    std::cout << "\tedge_vector = " << linearEdge->edgeVector(cubedSphere) << std::endl;
    std::cout << "\tsphLength = " << linearEdge->sphLength(cubedSphere) << std::endl;
    std::cout << "\teuclength = " << linearEdge->eucLength(cubedSphere) << std::endl;

    // EdgeSet
    EdgeSet<3> edges(linearEdgeFactory, SPHERICAL_SURFACE_GEOMETRY, 24);
    
    edges.insert(0, 1, 0, 1);
    edges.insert(1, 2, 0, 2);
    edges.divide(0, cubedSphere);
    edges.divide(1, cubedSphere);
    
    std::cout << edges.infoString(true) << std::endl;

#ifdef HAVE_VTK
	vtkSmartPointer<vtkCellArray> vcells = edges.toVtkCellArray();
	std::cout << "inserted " << vcells->GetNumberOfCells() << " cells to vtkCellArray." << std::endl;
#endif
    
return 0;
}

