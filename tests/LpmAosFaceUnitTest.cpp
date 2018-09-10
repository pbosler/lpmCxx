#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAosParticleFactory.hpp"
#include "LpmAosParticleSet.hpp"
#include "LpmAosEdgeFactory.hpp"
#include "LpmAosEdgeSet.hpp"
#include "LpmAosFace.hpp"
#include "LpmAosTriFace.hpp"
#include "LpmAosQuadFace.hpp"
#include "LpmAosQuadCubicFace.hpp"
#include "LpmAosFaceFactory.hpp"
#include "LpmAosFaceSet.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>
#ifdef HAVE_VTK

#endif

using namespace Lpm::Aos;
using Lpm::index_type;
using Lpm::scalar_type;
using Lpm::Logger;
using Lpm::OutputMessage;
using Lpm::GeometryType::SPHERICAL_SURFACE_GEOMETRY;
using Lpm::GeometryType::PLANAR_GEOMETRY;

typedef std::vector<index_type> ind_vec;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    std::stringstream ss;
    const std::string nullstr;
    {
        ss << "Test info: \n \t title: " << "Lpm Aos Face Unit Test" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify basic functionality of Lpm::Face class & its subclasses." << std::endl;
        ss << "\t 2. Verify basic functionality of Lpm::FaceSet class." << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
        ss.str(nullstr);
    }
    
    QuadFaceFactory<3> quadFactory;
    std::shared_ptr<FaceFactory<2>> triFactory(new TriFaceFactory<2>());
    
    const ind_vec ctr = {0};
    const ind_vec triverts = {0, 1, 2};
    const ind_vec quadverts = {0, 1, 2, 3};
    const ind_vec tedges = {10, 11, 12};
    const ind_vec quadedges = {10, 11, 12, 13};
    const index_type prt = 0;
    const scalar_type ar = 1.0;
    
    
    std::cout << "test start" << std::endl;
    std::cout << "Quad Face: " << std::endl;
    std::unique_ptr<Face<3>> qf = quadFactory.createFace(ctr, quadverts, quadedges, prt, ar);
    std::cout << qf->infoString();

    std::cout << "Tri face: " << std::endl;
    std::unique_ptr<Face<2>> tf = triFactory->createFace(ctr, triverts, tedges, prt, ar);
    std::cout << tf->infoString();

    
    std::shared_ptr<ParticleFactory<3>> sphereFactory(new SWEParticleFactory<3>());
    ParticleSet<3> cubedSphere(sphereFactory, 100);
    
    std::shared_ptr<ParticleFactory<2>> ppfac(new BasicParticleFactory<2>());
    ParticleSet<2> trihex(ppfac, 100);
    const scalar_type nullarea = 0.0;
	trihex.insert(Vec<2>(0.00000000000000000,0.00000000000000000), nullarea, true);
	trihex.insert(Vec<2>(0.50000000000000011,0.86602540378443860), nullarea, true);
	trihex.insert(Vec<2>(-0.49999999999999978,0.86602540378443871), nullarea, true);
	trihex.insert(Vec<2>(-1.00000000000000000,0.00000000000000000), nullarea, true);
	trihex.insert(Vec<2>(-0.50000000000000044,-0.86602540378443837), nullarea, true);
	trihex.insert(Vec<2>(0.49999999999999933,-0.86602540378443904), nullarea, true);
	trihex.insert(Vec<2>(1.00000000000000000,0.00000000000000000), nullarea, true);
	trihex.insert(Vec<2>(0.00000000000000011,0.57735026918962573));
	trihex.insert(Vec<2>(-0.49999999999999994,0.28867513459481292));
	trihex.insert(Vec<2>(-0.50000000000000011,-0.28867513459481281));
	trihex.insert(Vec<2>(-0.00000000000000037,-0.57735026918962584));
	trihex.insert(Vec<2>(0.49999999999999978,-0.28867513459481303));
	trihex.insert(Vec<2>(0.50000000000000000,0.28867513459481287));
	std::vector<std::vector<index_type>> edgeRecords(12);
	edgeRecords[0] = {0,1,0,5};
	edgeRecords[1] = {1,2,0,-1};
	edgeRecords[2] = {2,3,1,-1};
	edgeRecords[3] = {3,4,2,-1};
	edgeRecords[4] = {4,5,3,-1};
	edgeRecords[5] = {5,6,4,-1};
	edgeRecords[6] = {6,1,5,-1};
	edgeRecords[7] = {2,0,0,1};
	edgeRecords[8] = {3,0,1,2};
	edgeRecords[9] = {4,0,2,3};
	edgeRecords[10] = {0,5,4,3};
	edgeRecords[11] = {0,6,5,4};
	std::shared_ptr<EdgeFactory<2>> tefac(new LinearEdgeFactory<2>());
	EdgeSet<2> triedges(tefac, PLANAR_GEOMETRY, 100);
	for (int i=0; i<12; ++i) {
		triedges.insert(edgeRecords[i][0], edgeRecords[i][1], edgeRecords[i][2], edgeRecords[i][3]);
	}
	std::vector<std::vector<index_type>> tfv(6);
	tfv[0] = {0,1,2};
	tfv[1] = {2,3,0};
	tfv[2] = {4,0,3};
	tfv[3] = {0,4,5};
	tfv[4] = {5,6,0};
	tfv[5] = {1,0,6};
	std::vector<std::vector<index_type>> tfe(6);
	tfe[0] = {0,1,7};
	tfe[1] = {2,8,7};
	tfe[2] = {9,8,3};
	tfe[3] = {9,4,10};
	tfe[4] = {5,11,10};
	tfe[5] = {0,11,6};
	FaceSet<2> trifaces(triFactory, 100, PLANAR_GEOMETRY);
	const int root_index = -1;
	for (int i=0; i<6; ++i) {
		trifaces.insert(ind_vec(1, i+7), tfv[i], tfe[i], root_index);
	}
	trifaces.setArea(trihex);
	std::cout << trifaces.infoString() << std::endl;
	std::cout << triedges.infoString() << std::endl;
	std::cout << trihex.infoString(true) << std::endl;
	std::cout << "calling TriFace::divide " << std::endl;
	trifaces.divide(0, trihex, triedges);
	std::cout << trifaces.infoString(true) << std::endl;
	std::cout << triedges.infoString(true) << std::endl;
	std::cout << trihex.infoString(true) << std::endl;
	
	
    
    
 //    cubedSphere.insert(Vec<3>(0.57735026918962584, -0.57735026918962584,  0.57735026918962584));
//     cubedSphere.insert(Vec<3>(0.57735026918962584, -0.57735026918962584, -0.57735026918962584));
//     cubedSphere.insert(Vec<3>(0.57735026918962584,  0.57735026918962584, -0.57735026918962584));
//     cubedSphere.insert(Vec<3>(0.57735026918962584,  0.57735026918962584,  0.57735026918962584));
//     cubedSphere.insert(Vec<3>(-0.57735026918962584,  0.57735026918962584, -0.57735026918962584));
//     cubedSphere.insert(Vec<3>(-0.57735026918962584,  0.57735026918962584,  0.57735026918962584));
//     cubedSphere.insert(Vec<3>(-0.57735026918962584, -0.57735026918962584, -0.57735026918962584));
//     cubedSphere.insert(Vec<3>(-0.57735026918962584, -0.57735026918962584,  0.57735026918962584));
//     cubedSphere.insert(Vec<3>(1.0000000000000000,   0.0000000000000000,   0.0000000000000000), 2.0943951023931948);
//     cubedSphere.insert(Vec<3>(0.0000000000000000,   1.0000000000000000,   0.0000000000000000), 2.0943951023931948);
//     cubedSphere.insert(Vec<3>(-1.000000000000000,   0.0000000000000000,   0.0000000000000000), 2.0943951023931948);
//     cubedSphere.insert(Vec<3>(0.0000000000000000,  -1.0000000000000000,   0.0000000000000000), 2.0943951023931948);
//     cubedSphere.insert(Vec<3>(0.0000000000000000,   0.0000000000000000,   1.0000000000000000), 2.0943951023931948);
//     cubedSphere.insert(Vec<3>(0.0000000000000000,   0.0000000000000000,  -1.0000000000000000), 2.0943951023931948);
// 
//     std::cout << cubedSphere.infoString() << std::endl;
//     
//     
//     edgeRecords[0] = {0,1,0,3};
//     edgeRecords[1] = {1,2,0,5};
//     edgeRecords[2] = {2,3,0,1};
// 	edgeRecords[3] = {3,0,0,4};
// 	edgeRecords[4] = {2,4,1,5};
// 	edgeRecords[5] = {4,5,1,2};
// 	edgeRecords[6] = {5,3,1,4};
// 	edgeRecords[7] = {4,6,2,5};
// 	edgeRecords[8] = {6,7,2,3};
// 	edgeRecords[9] = {7,5,2,4};
// 	edgeRecords[10] = {6,1,3,5};
// 	edgeRecords[11] = {0,7,3,4};
// 	
// 	std::shared_ptr<EdgeFactory<3>> edgeFac(new LinearEdgeFactory<3>());
// 	EdgeSet<3> edges(edgeFac, SPHERICAL_SURFACE_GEOMETRY, 100);
// 	for (int i=0; i<12; ++i) {
// 		edges.insert(edgeRecords[i][0], edgeRecords[i][1], edgeRecords[i][2], edgeRecords[i][3]);
// 	}
// 	std::cout << edges.infoString();
// 	
// 	
//     std::vector<std::vector<index_type>> fv(6, std::vector<index_type>(4));
//     fv[0] = {1, 2, 3, 4};
//     fv[1] = {4, 3, 5, 6};
//     fv[2] = {6, 5, 7, 8};
//     fv[3] = {8, 7, 2, 1};
//     fv[4] = {8, 1, 4, 6};
//     fv[5] = {2, 7, 5, 3};
//     std::vector<std::vector<index_type>> fe(6, std::vector<index_type>(4));
//     fe[0] = {1, 2, 3, 4};
//     fe[1] = {3, 5, 6, 7};
//     fe[2] = {6, 8, 9, 10};
//     fe[3] = {9, 11, 1, 12};
//     fe[4] = {12, 4, 7, 10};
//     fe[5] = {11, 8, 5, 2};
//     for (int i=0; i<6; ++i) {
//         for (int j=0; j<4; ++j) {
//             fv[i][j] -= 1;
//             fe[i][j] -= 1;
//         }
//     }
//     std::shared_ptr<FaceFactory<3>> csFaceFac = std::shared_ptr<FaceFactory<3>>(new QuadFaceFactory<3>());
//     FaceSet<3> cs(csFaceFac, 100, SPHERICAL_SURFACE_GEOMETRY);
// 
//     for (int i=0; i<6; ++i) {
//         cs.insert(ind_vec(1, i+8), fv[i], fe[i], root_index);
//     }
//     std::cout << cs.infoString() << std::endl;
//     cs.setArea(cubedSphere);
//     std::cout << "calling divide" << std::endl;
//     cs.divide(0, cubedSphere, edges);
// 
// //     std::cout << "edges:" << std::endl;
// //     std::cout << edges.infoString(true) << std::endl;
// //     std::cout << cs.infoString(true) << std::endl;
// //     cs.setArea(cubedSphere);
// //     std::cout << "returned from set area" << std::endl;
//     std::cout << cs.infoString() << std::endl;
//     
//     
// #ifdef HAVE_VTK
// 	std::cout << "Testing VTK" << std::endl;
// 	vtkSmartPointer<vtkCellArray> vcells = cs.toVtkCellArray();
// 	std::cout << "inserted " << vcells->GetNumberOfCells() << " cells to vtkCellArray." << std::endl;
// #endif    
    
return 0;
}

