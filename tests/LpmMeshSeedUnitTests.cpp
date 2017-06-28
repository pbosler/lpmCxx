#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmMeshSeed.h"
#include "LpmTriFaces.h"
#include "LpmQuadFaces.h"
#include "LpmEuclideanCoords.h"
#include "LpmSphericalCoords.h"
#include <memory>
#include <iostream>

using namespace Lpm;

int main(int argc, char* argv[]) {
    {
    std::cout << "**** TriHexSeed ****" << std::endl;
    std::shared_ptr<EuclideanCoords> ec(new EuclideanCoords(20));
    std::shared_ptr<Edges> edges(new Edges(50, ec));
    std::shared_ptr<TriFaces> faces(new TriFaces(12, edges, ec));
    
    TriHexSeed mSeed;
    mSeed.initMeshFromSeed(ec, edges, faces);
    std::cout << mSeed.MeshSeed::infoString() << std::endl;
    }
    {
    std::cout << "**** IcosTriSphereSeed ****" << std::endl;
    std::shared_ptr<SphericalCoords> ec(new SphericalCoords(12));
    std::shared_ptr<Edges> edges(new Edges(30, ec));
    std::shared_ptr<TriFaces> faces(new TriFaces(20, edges, ec));
    
    IcosTriSphereSeed mSeed;
    mSeed.initMeshFromSeed(ec, edges, faces);
    std::cout << mSeed.MeshSeed::infoString() << std::endl;
    }
    {
    std::cout << "**** QuadRectSeed ****" << std::endl;
    std::shared_ptr<EuclideanCoords> ec(new EuclideanCoords(20));
    std::shared_ptr<Edges> edges(new Edges(50, ec));
    std::shared_ptr<TriFaces> faces(new TriFaces(12, edges, ec));
    
    QuadRectSeed mSeed(1.0);
    mSeed.initMeshFromSeed(ec, edges, faces);
    std::cout << mSeed.MeshSeed::infoString() << std::endl;
    }
    {
    std::cout << "**** CubedSphereSeed ****" << std::endl;
    std::shared_ptr<SphericalCoords> ec(new SphericalCoords(20));
    std::shared_ptr<Edges> edges(new Edges(50, ec));
    std::shared_ptr<TriFaces> faces(new TriFaces(12, edges, ec));
    
    CubedSphereSeed mSeed;
    mSeed.initMeshFromSeed(ec, edges, faces);
    std::cout << mSeed.MeshSeed::infoString() << std::endl;
    }    
return 0;
}