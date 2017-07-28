#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOutputMessage.h"
#include "LpmLogger.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"
#include "LpmBox3d.h"
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
        ss << "\t 1. Verify Box3d class functions" << std::endl;
        ss << "\t 2. Verify basic Treenode class functions" << std::endl;
        ss << "\t 3. Generate 3d octree for Cartesian points" << std::endl;
        ss << "\t 4. Generate 3d octree for spherical points" << std::endl;
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
            
    }
    {
        const int nMax = 8000;
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
        std::shared_ptr<Tree> tree(new Tree(ec, maxAspectRatio));
        std::cout << "tree info: " << tree->infoString();

        const int nCoordsPerNode = 10;
        std::cout << "calling generateTree..." << std::endl;
        tree->buildTree(nCoordsPerNode);
        std::cout << "returned from generateTree:" << std::endl;
        std::cout << "\t nNodes = " << tree->nNodes() << std::endl;
        std::cout << "\t treeDepth = " << tree->depth() << std::endl;
        

        
        const std::string fname("octreeUnitTest.vtk");
        std::stringstream ss;
        ss << "nCoords = " << nMax << ", nCoordsPerNode = " << nCoordsPerNode;
        tree->writeToVtk(fname, ss.str());
    }
    
    {
        const int nMax = 25000;
        std::shared_ptr<SphericalCoords> sc(new SphericalCoords(nMax));
        GeometryType geom = sc->geometry();
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
        
        
        sc->initRandom();
        
        std::ofstream cfile("octreeSphereCoords.txt");
        sc->writeCoordsCSV(cfile);
        cfile.close();
        
        const scalar_type maxAspectRatio = 2.0;
        
        std::shared_ptr<Tree> tree(new Tree(sc, maxAspectRatio));
        std::cout << tree->infoString();
        
        const int nCoordsPerNode = 20;
        
        tree->buildTree(nCoordsPerNode);
        std::cout << "returned from generateTree:" << std::endl;
        std::cout << "\t nNodes = " << tree->nNodes() << std::endl;
        std::cout << "\t nRecursiveNodes " << tree->recursiveNNodes() << std::endl;
        std::cout << "\t treeDepth = " << tree->depth() << std::endl;
        
        std::cout << tree->infoString();
        
        const std::string fname("sphereOctreeUnitTest.vtk");
        std::stringstream ss;
        ss << "nCoords = " << nMax << ", nCoordsPerNode = " << nCoordsPerNode;
        tree->writeToVtk(fname, ss.str());
    }

return 0;
}

