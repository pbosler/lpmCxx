#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmAosFace.hpp"
#include "LpmAosFaceFactory.hpp"
#include "LpmAosFaceSet.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>

using namespace Lpm;

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
    TriFaceFactory<2> triFactory;
    
    const ind_vec ctr = {0};
    const ind_vec triverts = {0, 1, 2};
    const ind_vec quadverts = {0, 1, 2, 3};
    const ind_vec triedges = {10, 11, 12};
    const ind_vec quadedges = {10, 11, 12, 13};
    const index_type prt = 0;
    const scalar_type ar = 1.0;
    
    
    std::cout << "test start" << std::endl;
    std::cout << "Quad Face: " << std::endl;
    std::unique_ptr<Face<3>> qf = quadFactory.createFace(ctr, quadverts, quadedges, prt, ar);
    std::cout << qf->infoString();

    std::cout << "Tri face: " << std::endl;
    std::unique_ptr<Face<2>> tf = triFactory.createFace(ctr, triverts, triedges, prt, ar);
    std::cout << tf->infoString();
    
    std::shared_ptr<ParticleFactory<3>> sphereFactory(new SWEParticleFactory<3>());
    ParticleSet<3> cubedSphere(sphereFactory, 15);
    
    cubedSphere.insert(Vec<3>(0.57735026918962584, -0.57735026918962584,  0.57735026918962584));
    cubedSphere.insert(Vec<3>(0.57735026918962584, -0.57735026918962584, -0.57735026918962584));
    cubedSphere.insert(Vec<3>(0.57735026918962584,  0.57735026918962584, -0.57735026918962584));
    cubedSphere.insert(Vec<3>(0.57735026918962584,  0.57735026918962584,  0.57735026918962584));
    cubedSphere.insert(Vec<3>(-0.57735026918962584,  0.57735026918962584, -0.57735026918962584));
    cubedSphere.insert(Vec<3>(-0.57735026918962584,  0.57735026918962584,  0.57735026918962584));
    cubedSphere.insert(Vec<3>(-0.57735026918962584, -0.57735026918962584, -0.57735026918962584));
    cubedSphere.insert(Vec<3>(-0.57735026918962584, -0.57735026918962584,  0.57735026918962584));
    cubedSphere.insert(Vec<3>(1.0000000000000000,   0.0000000000000000,   0.0000000000000000), 2.0943951023931948);
    cubedSphere.insert(Vec<3>(0.0000000000000000,   1.0000000000000000,   0.0000000000000000), 2.0943951023931948);
    cubedSphere.insert(Vec<3>(-1.000000000000000,   0.0000000000000000,   0.0000000000000000), 2.0943951023931948);
    cubedSphere.insert(Vec<3>(0.0000000000000000,  -1.0000000000000000,   0.0000000000000000), 2.0943951023931948);
    cubedSphere.insert(Vec<3>(0.0000000000000000,   0.0000000000000000,   1.0000000000000000), 2.0943951023931948);
    cubedSphere.insert(Vec<3>(0.0000000000000000,   0.0000000000000000,  -1.0000000000000000), 2.0943951023931948);

    std::cout << cubedSphere.infoString() << std::endl;

    std::vector<std::vector<index_type>> fv(6, std::vector<index_type>(4));
    fv[0] = {1, 2, 3, 4};
    fv[1] = {4, 3, 5, 6};
    fv[2] = {6, 5, 7, 8};
    fv[3] = {8, 7, 2, 1};
    fv[4] = {8, 1, 4, 6};
    fv[5] = {2, 7, 5, 3};
    std::vector<std::vector<index_type>> fe(6, std::vector<index_type>(4));
    fe[0] = {1, 2, 3, 4};
    fe[1] = {3, 5, 6, 7};
    fe[2] = {6, 8, 9, 10};
    fe[3] = {9, 11, 1, 12};
    fe[4] = {12, 4, 7, 10};
    fe[5] = {11, 8, 5, 2};
    for (int i=0; i<6; ++i) {
        for (int j=0; j<4; ++j) {
            fv[i][j] -= 1;
            fe[i][j] -= 1;
        }
    }
    std::shared_ptr<FaceFactory<3>> csFaceFac = std::shared_ptr<FaceFactory<3>>(new QuadFaceFactory<3>());
    FaceSet<3> cs(csFaceFac, 6, SPHERICAL_SURFACE_GEOMETRY);
    for (int i=0; i<6; ++i) {
        cs.insert(ind_vec(1, i+8), fv[i], fe[i], -1);
    }
    
    std::cout << cs.infoString() << std::endl;
return 0;
}

