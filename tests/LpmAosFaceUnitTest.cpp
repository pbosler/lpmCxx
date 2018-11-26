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
#include "LpmGll.hpp"
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
    {
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
	
	}
	{ // Cubic section
		
		std::cout << std::endl << "vvvvvvvvvv CUBIC SECTION vvvvvvvvvv" << std::endl;
	
		std::shared_ptr<ParticleFactory<2>> cpfac(new SWEParticleFactory<2>());
		std::shared_ptr<EdgeFactory<2>> cefac(new CubicEdgeFactory<2>());
		std::shared_ptr<FaceFactory<2>> cffac(new QuadCubicFaceFactory<2>());
		
		Lpm::CubicGLL<2> gll;
		
		ParticleSet<2> cps(cpfac, 100);
		const scalar_type map_jac = 0.25;
		std::vector<Vec<2>> corners(4);
		corners[0] = Vec<2>(-1,1);
		corners[1] = Vec<2>(-1,0);
		corners[2] = Vec<2>(0,0);
		corners[3] = Vec<2>(0,1);
		for (short i=0; i<12; ++i)
			cps.insert(gll.bilinearMap(corners, gll.quad16edgeqp(i)), map_jac*gll.quad16edgeqw(i), true);
		for (short i=0; i<4; ++i)
			cps.insert(gll.bilinearMap(corners, gll.quad16centerqp(i)), map_jac*gll.quad16centerqw(i), false);
		corners[0] = Vec<2>(-1,0);
		corners[1] = Vec<2>(-1,-1);
		corners[2] = Vec<2>(0,-1);
		corners[3] = Vec<2>(0,0);
		for (short i=1; i<=8; ++i)
			cps.insert(gll.bilinearMap(corners, gll.quad16edgeqp(i)), map_jac*gll.quad16edgeqw(i), true);
		for (short i=0; i<4; ++i)
			cps.insert(gll.bilinearMap(corners, gll.quad16centerqp(i)), map_jac*gll.quad16centerqw(i), false);
		corners[0] = Vec<2>(0,0);
		corners[1] = Vec<2>(0,-1);
		corners[2] = Vec<2>(1,-1);
		corners[3] = Vec<2>(1,0);
		for (short i=4; i<=11; ++i)
			cps.insert(gll.bilinearMap(corners, gll.quad16edgeqp(i)), map_jac*gll.quad16edgeqw(i), true);
		for (short i=0; i<4; ++i)
			cps.insert(gll.bilinearMap(corners, gll.quad16centerqp(i)), map_jac*gll.quad16centerqw(i), false);
		corners[0] = Vec<2>(0,1);
		corners[1] = Vec<2>(0,0);
		corners[2] = Vec<2>(1,0);
		corners[3] = Vec<2>(1,1);
		for (short i=7; i<=11; ++i)
			cps.insert(gll.bilinearMap(corners, gll.quad16edgeqp(i)), map_jac*gll.quad16edgeqw(i), true);
		for (short i=0; i<4; ++i)
			cps.insert(gll.bilinearMap(corners, gll.quad16centerqp(i)), map_jac*gll.quad16centerqw(i), false);
		
		
// 		std::cout << cps.infoString(true);
// 		for (int i=0; i<49; ++i)
// 			std::cout << i << ": " << cps.weight(i) << std::endl;
				
		EdgeSet<2> ces(cefac, PLANAR_GEOMETRY, 100);
		std::vector<std::vector<index_type>> edge_ints(12,std::vector<index_type>(2,-1));
		edge_ints[0] = {1,2};
		edge_ints[1] = {16,17};
		edge_ints[2] = {19,20};
		edge_ints[3] = {28,29};
		edge_ints[4] = {31,32};
		edge_ints[5] = {40,41};
		edge_ints[6] = {43,44};
		edge_ints[7] = {10,11};
		edge_ints[8] = {4,5};
		edge_ints[9] = {35,35};
		edge_ints[10] = {22,23};
		edge_ints[11] = {7,8};
		std::vector<index_type> edgeOrigs = {0,3,18,21,30,33,42,9,3,6,21,6};
		std::vector<index_type> edgeDests ={3,18,21,30,33,42,9,0,6,33,6,9};
		std::vector<index_type> edgeLefts = {0,1,1,2,2,3,3,0,0,3,1,0};
		std::vector<index_type> edgeRights = {-1,-1,-1,-1,-1,-1,-1,-1,1,2,2,3};
		for (int i=0; i<12; ++i) 
			ces.insert(edgeOrigs[i], edgeDests[i], edgeLefts[i], edgeRights[i], edge_ints[i]);
		
		FaceSet<2> cfs(cffac, 100, PLANAR_GEOMETRY);
		std::vector<std::vector<index_type>> fverts(4);
		fverts[0] = {0,1,2,3,4,5,6,7,8,9,10,11};
		fverts[1] = {3,16,17,18,19,20,21,22,23,6,5,4};
		fverts[2] = {6,23,22,21,28,29,30,31,32,33,34,35};
		fverts[3] = {9,8,7,6,35,34,33,40,41,42,43,44};
		
		std::vector<std::vector<index_type>> fedges(4);
		fedges[0] = {0,8,11,7};
		fedges[1] = {1,2,10,8};
		fedges[2] = {10,3,4,9};
		fedges[3] = {11,9,5,6};
		
		std::vector<std::vector<index_type>> fintrs(4);
		fintrs[0] = {12,13,14,15};
		fintrs[1] = {24,25,26,27};
		fintrs[2] = {36,37,38,39};
		fintrs[3] = {45,46,47,48};
		
		const index_type parent_index = -1;
		const scalar_type root_area = 1.0;
		for (int i=0; i<4; ++i) {
			cfs.insert(fintrs[i], fverts[i], fedges[i], parent_index, root_area);
		}
		
		std::cout << "CUBIC PLANE INFO" << std::endl;
		std::cout << cps.infoString();
		std::cout << ces.infoString();
		std::cout << cfs.infoString();
		
		std::cout << "CUBIC PLANE INFO  -- after divide face 0" << std::endl;
		cfs.divide(0, cps, ces);
		std::cout << cps.infoString();
		std::cout << ces.infoString();
		std::cout << cfs.infoString();
		
		std::cout << "RESET AREA" << std::endl;
		cfs.setArea(cps);
		std::cout << cfs.infoString();
	}

    
return 0;
}

