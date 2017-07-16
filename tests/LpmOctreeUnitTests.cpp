#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"
#include "LpmOctree.h"
#include "LpmCoords.h"
#include "LpmXyzVector.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <numeric>

using namespace Lpm;

int main (int argc, char* argv[]) {

    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    {
        std::stringstream ss;
        ss << "Test info: \n \t title: " << "Octree unit tests" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Generate 3d octree for Cartesian points" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
    }
    {   // box 3d unit tests
        
        Box3d unitBox(-1,1,-1,1,-1,1);
        std::cout << "Unit Box : " << unitBox.infoString();
        std::cout << "\tvolume = " << unitBox.volume() << std::endl;
        std::cout << "\tcentroid = " << unitBox.centroid() << std::endl;
        std::cout << "\tcontains origin ? : " << (unitBox.containsPoint(XyzVector(0.0, 0.0, 0.0)) ? "true" : "false") << std::endl;
        std::cout << "\tlongestEdge = " << unitBox.longestEdge() << std::endl;
        std::cout << "\tshortestEdge = " << unitBox.shortestEdge() << std::endl;
        std::cout << "\taspectRatio = " << unitBox.aspectRatio() << std::endl;
        
        std::vector<Box3d> kids = unitBox.bisectAll();
        for (int i = 0; i < 8; ++i) {
            std::cout << "child " << i << ": " << kids[i].infoString();
        }
        
        Box3d boxCopy(kids[1]);
        std::cout << "BoxCopy " << boxCopy.infoString();
        Box3d boxAssign = boxCopy;
        std::cout << "BoxAssign " << boxAssign.infoString();
        
        std::shared_ptr<Treenode> nullnode(new Treenode());
        std::cout << "nullnode info: ";
        std::cout << nullnode->infoString();
        std::vector<index_type> inds = {0,1,2,3};
        std::vector<std::shared_ptr<Treenode>> pvec;
        nullnode->children.push_back(std::shared_ptr<Treenode>(new Treenode(kids[1], nullnode, inds, 1.0)));
        std::cout << "nullnode child 0: " ;
        std::cout << nullnode->children[0]->infoString();
    
    }
    {
        const int nMax = 1000;
        const scalar_type domainRadius = 2.0;
        
        std::shared_ptr<EuclideanCoords> ec(new EuclideanCoords(nMax, CARTESIAN_3D_GEOMETRY));
        GeometryType geom = ec->geometry();
        std::cout << "geometry: " << geom << std::endl;
        switch (geom) {
            case PLANAR_GEOMETRY : {
                std::cout << "planar geometry";
                break;
            }
            case SPHERICAL_SURFACE_GEOMETRY : {
                std::cout << "spherical surface geometry";
                break;
            }
            case CARTESIAN_3D_GEOMETRY : {
                std::cout << "Cartesian 3D geometry";
            }
        }
        std::cout << std::endl;
        
        ec->initRandom(false, domainRadius);
        std::ofstream cfile("octreeCoords.txt");
        ec->writeCoordsCSV(cfile);
        cfile.close();
        
        const scalar_type maxAspectRatio = 1.5;
        std::shared_ptr<Treenode> tree(new Treenode(ec, maxAspectRatio));
        std::cout << "tree info: " << tree->infoString();
        
        std::cout << "Root node: " << std::endl;
        std::cout << "\t" << tree->box.infoString();
        
        
        const int nCoordsPerNode = 10;
        std::cout << "calling generateTree..." << std::endl;
        generateTree(tree, ec, nCoordsPerNode);
        std::cout << "returned from generateTree; nNodes = " << nTreenodes(tree) << std::endl;
        

        
        const std::string fname("octreeUnitTest.vtk");
        std::stringstream ss;
        ss << "nCoords = " << nMax << ", nCoordsPerNode = " << nCoordsPerNode;
        writeTreeToVtk(fname, ss.str(), tree);
    }
    
    {
        const int nMax = 5000;
        std::shared_ptr<SphericalCoords> sc(new SphericalCoords(nMax));
        
        sc->initRandom();
        
        std::ofstream cfile("octreeSphereCoords.txt");
        sc->writeCoordsCSV(cfile);
        cfile.close();
        
        const scalar_type maxAspectRatio = 2.0;
        
        std::shared_ptr<Treenode> tree(new Treenode(sc, maxAspectRatio));
        
        const int nCoordsPerNode = 20;
        
        generateTree(tree, sc, nCoordsPerNode);
        
        const std::string fname("sphereOctreeUnitTest.vtk");
        std::stringstream ss;
        ss << "nCoords = " << nMax << ", nCoordsPerNode = " << nCoordsPerNode;
        writeTreeToVtk(fname, ss.str(), tree);
    }

return 0;
}

