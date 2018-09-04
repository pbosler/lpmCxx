#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAosPolyMesh2d.hpp"
// #include "LpmAosMeshSeedFactory.hpp"
#include "LpmAosParticleFactory.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>

#ifdef HAVE_VTK
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#endif

using namespace Lpm::Aos;
using Lpm::scalar_type;
using Lpm::index_type;
using Lpm::Logger;
using Lpm::OutputMessage;
using Lpm::geometryString;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "LpmAosPolyMesh2d Unit tests" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify initial meshes from seed." << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }
    
    std::shared_ptr<MeshSeed<2>> triPlaneSeed = std::shared_ptr<MeshSeed<2>>(new TriHexSeed());
    std::shared_ptr<MeshSeed<2>> quadPlaneSeed = std::shared_ptr<MeshSeed<2>>(new QuadRectSeed());
    std::shared_ptr<MeshSeed<3>> triSphereSeed = std::shared_ptr<MeshSeed<3>>(new IcosTriSphereSeed());
    std::shared_ptr<MeshSeed<3>> quadSphereSeed = std::shared_ptr<MeshSeed<3>>(new CubedSphereSeed());
    
    std::shared_ptr<BasicParticleFactory<2>> pfac_plane(new BasicParticleFactory<2>());
    std::shared_ptr<SWEParticleFactory<3>> pfac_sphere(new SWEParticleFactory<3>());
    
    std::shared_ptr<LinearEdgeFactory<2>> efac_plane(new LinearEdgeFactory<2>());
    std::shared_ptr<LinearEdgeFactory<3>> efac_sphere(new LinearEdgeFactory<3>());
    
    std::shared_ptr<TriFaceFactory<2>> trifac_plane(new TriFaceFactory<2>());
    std::shared_ptr<TriFaceFactory<3>> trifac_sphere(new TriFaceFactory<3>());
    std::shared_ptr<QuadFaceFactory<2>> quadfac_plane(new QuadFaceFactory<2>());
    std::shared_ptr<QuadFaceFactory<3>> quadfac_sphere(new QuadFaceFactory<3>());

    const int initnest = 1;
    const int maxnest = initnest;
    index_type nMaxTriPlaneParticles;
    index_type nMaxTriPlaneEdges;
    index_type nMaxTriPlaneFaces;
    index_type nMaxQuadPlaneParticles;
    index_type nMaxQuadPlaneEdges;
    index_type nMaxQuadPlaneFaces;
    index_type nMaxTriSphereParticles;
    index_type nMaxTriSphereEdges;
    index_type nMaxTriSphereFaces;
    index_type nMaxQuadSphereParticles;
    index_type nMaxQuadSphereEdges;
    index_type nMaxQuadSphereFaces;
    const int amrLimit = 0;
    const scalar_type radius = 1.0;

    triPlaneSeed->determineMaxAllocations(nMaxTriPlaneParticles, nMaxTriPlaneEdges, nMaxTriPlaneFaces, maxnest);
    quadPlaneSeed->determineMaxAllocations(nMaxQuadPlaneParticles, nMaxQuadPlaneEdges, nMaxQuadPlaneFaces, maxnest);
    triSphereSeed->determineMaxAllocations(nMaxTriSphereParticles, nMaxTriSphereEdges, nMaxTriSphereFaces, maxnest);
    quadSphereSeed->determineMaxAllocations(nMaxQuadSphereParticles, nMaxQuadSphereEdges, nMaxQuadSphereFaces, maxnest);
    
    PolyMesh2d<2> triplane(triPlaneSeed, pfac_plane, efac_plane, trifac_plane, initnest, maxnest, amrLimit, radius);
    triplane.initFromSeedStaggeredFacesAndVerts();
    
    PolyMesh2d<2> quadplane(quadPlaneSeed, pfac_plane, efac_plane, quadfac_plane, initnest, maxnest, amrLimit, radius);
    quadplane.initFromSeedStaggeredFacesAndVerts();
    
    PolyMesh2d<3> trisphere(triSphereSeed, pfac_sphere, efac_sphere, trifac_sphere, initnest, maxnest, amrLimit, radius);
    trisphere.initFromSeedStaggeredFacesAndVerts();
    
    PolyMesh2d<3> quadsphere(quadSphereSeed, pfac_sphere, efac_sphere, quadfac_sphere, initnest, maxnest, amrLimit, radius);
    quadsphere.initFromSeedStaggeredFacesAndVerts();
    
    std::cout << "all meshes with initnest = " << initnest << " created." << std::endl;
    
    std::cout << std::endl << std::endl;
    std::cout << triplane.infoString();
    std::cout << std::endl << std::endl;
    
    std::cout << quadplane.infoString();
    std::cout << std::endl << std::endl;
    
    std::cout << trisphere.infoString();
    std::cout << std::endl <<std::endl;
    
    std::cout << quadsphere.infoString();
    std::cout << std::endl << std::endl;
    
#ifdef HAVE_VTK
    {
    const std::string froot = "tmp/polymeshtest_";
    std::ostringstream ss;
    ss << froot << "planetri_" << initnest << ".vtk";
    vtkSmartPointer<vtkPolyData> pd = triplane.toVtkPolyData();
    std::cout << "returned from vtkPolyData conversion." << std::endl;
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(pd);
    writer->SetFileName(ss.str().c_str());
    writer->Write();
    ss.str(std::string());
    
    ss << froot << "planequad_" << initnest << ".vtk";
    pd = quadplane.toVtkPolyData();
    writer->SetInputData(pd);
    writer->SetFileName(ss.str().c_str());
    writer->Write();
    ss.str(std::string());
    
    ss << froot << "trisphere_" << initnest << ".vtk";
    pd = trisphere.toVtkPolyData();
    writer->SetInputData(pd);
    writer->SetFileName(ss.str().c_str());
    writer->Write();
    ss.str(std::string());
    
    ss << froot << "quadsphere" << initnest << ".vtk";
    pd = quadsphere.toVtkPolyData();
    writer->SetInputData(pd);
    writer->SetFileName(ss.str().c_str());
    writer->Write();
    ss.str(std::string());
    }
#endif    
    
return 0;
}

