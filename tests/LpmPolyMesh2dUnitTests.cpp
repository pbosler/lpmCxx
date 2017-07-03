#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmPolyMesh2d.h"
#include "LpmEdges.h"
#include "LpmFaces.h"
#include "LpmTriFaces.h"
#include "LpmQuadFaces.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"
#include "LpmPolyMesh2d.h"
#include "LpmVtkFileIO.h"
#include "LpmXyzVector.h"
#include <memory>
#include <iostream>
#include <sstream>

using namespace Lpm;

int main (int argc, char* argv[]) {
    std::unique_ptr<Logger> log(new Logger(OutputMessage::debugPriority));
    {
        std::stringstream ss;
        ss << "Test info: \n \t title: " << "PolyMesh2d Unit Tests" << std::endl;
        ss << "\t objectives: " << std::endl;
        ss << "\t 1. Verify PolyMesh2d constructors, .vtk output" << std::endl;
        ss << "\t 2. Verify basic topology operations (*TO DO*)" << std::endl;
        OutputMessage introMsg(ss.str(), OutputMessage::tracePriority, "main");
        log->logMessage(introMsg);
    }
    const int maxRecursion = 1;
    const scalar_type dradius = 3.0;
    
    {
        std::cout << "**** Planar triangle mesh ****" << std::endl;
        TriHexSeed mSeed;
        const bool isLagrangian = false;
        PolyMesh2d mesh(mSeed, maxRecursion, isLagrangian, dradius);
        
        std::stringstream ss;
        ss << mSeed.idString() << maxRecursion << ".vtk";
        mesh.writeToVTKFile(ss.str(), "unitTest");
        
        const XyzVector queryPt(0.2, 0.05);
        const index_type face_index = mesh.locateFaceContainingPoint(queryPt);
        const std::vector<XyzVector> enclosingVertices = mesh.getCoordVecs(mesh.faceVertices(face_index));
        std::cout << "found query point " << queryPt << " in face " << face_index << std::endl;
        std::cout << "\tface " << face_index << " has vertex coordinates:" << std::endl;
        for (index_type i = 0; i < enclosingVertices.size(); ++i)
            std::cout << "\t" << enclosingVertices[i] << std::endl;
    }
    {
        std::cout << "**** Spherical triangle mesh ****" << std::endl;
    
        IcosTriSphereSeed mSeed;
        const bool isLagrangian = true;
        PolyMesh2d mesh(mSeed, maxRecursion, isLagrangian, dradius);
    
        std::stringstream ss;
        ss << mSeed.idString() << maxRecursion << ".vtk";
        mesh.writeToVTKFile(ss.str(), "unitTest");
        
        const scalar_type lam0 = 0.2;
        const scalar_type the0 = 0.05;
        scalar_type xx, yy, zz;
        llToXyz(xx, yy, zz, lam0, the0);
        const XyzVector queryPt(xx, yy, zz);
        const index_type face_index = mesh.locateFaceContainingPoint(queryPt);
        const std::vector<XyzVector> enclosingVertices = mesh.getCoordVecs(mesh.faceVertices(face_index));
        std::cout << "found query point " << queryPt << " in face " << face_index << std::endl;
        std::cout << "\tface " << face_index << " has vertex coordinates:" << std::endl;
        for (index_type i = 0; i < enclosingVertices.size(); ++i)
            std::cout << "\t" << enclosingVertices[i] << std::endl;
    }
    {
        std::cout << "**** Planar quadrilateral mesh ****" << std::endl;
    
        QuadRectSeed mSeed;
        const bool isLagrangian = true;
        PolyMesh2d mesh(mSeed, maxRecursion, isLagrangian, dradius);
        
        std::stringstream ss;
        ss << mSeed.idString() << maxRecursion << ".vtk";
        mesh.writeToVTKFile(ss.str(), "unitTest");
        
        const XyzVector queryPt(0.2, 0.05);
        const index_type face_index = mesh.locateFaceContainingPoint(queryPt);
        const std::vector<XyzVector> enclosingVertices = mesh.getCoordVecs(mesh.faceVertices(face_index));
        std::cout << "found query point " << queryPt << " in face " << face_index << std::endl;
        std::cout << "\tface " << face_index << " has vertex coordinates:" << std::endl;
        for (index_type i = 0; i < enclosingVertices.size(); ++i)
            std::cout << "\t" << enclosingVertices[i] << std::endl;
    }
    {
        std::cout << "**** Spherical quadrilateral mesh ****" << std::endl;
    
        CubedSphereSeed mSeed;
        const bool isLagrangian = false;
        PolyMesh2d mesh(mSeed, maxRecursion, isLagrangian, dradius);
        
        std::stringstream ss;
        ss << mSeed.idString() << maxRecursion << ".vtk";
        mesh.writeToVTKFile(ss.str(), "unitTest");
        
        const scalar_type lam0 = 0.2;
        const scalar_type the0 = 0.05;
        scalar_type xx, yy, zz;
        llToXyz(xx, yy, zz, lam0, the0);
        const XyzVector queryPt(xx, yy, zz);
        const index_type face_index = mesh.locateFaceContainingPoint(queryPt);
        const std::vector<XyzVector> enclosingVertices = mesh.getCoordVecs(mesh.faceVertices(face_index));
        std::cout << "found query point " << queryPt << " in face " << face_index << std::endl;
        std::cout << "\tface " << face_index << " has vertex coordinates:" << std::endl;
        for (index_type i = 0; i < enclosingVertices.size(); ++i)
            std::cout << "\t" << enclosingVertices[i] << std::endl;
    }  

return 0;
}

