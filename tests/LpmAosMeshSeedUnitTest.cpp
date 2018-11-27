#include "LpmAosMeshSeed.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

using namespace Lpm::Aos;
using Lpm::index_type;

int main(int argc, char* argv[]) {

    index_type nv;
    index_type ne;
    index_type nf;
    
    const index_type maxRec = 9;

    TriHexSeed triPlane;
    std::cout << "id string = " << triPlane.idString() << std::endl;
    std::cout << "returned." << std::endl;
    std::cout << triPlane.infoString();
    
    QuadRectSeed quadPlane;
    std::cout << quadPlane.infoString();
    
    IcosTriSphereSeed icosTriSphere;\
    std::cout << icosTriSphere.infoString();
    
    CubedSphereSeed cubedSphere;
    std::cout << cubedSphere.infoString();
    
    QuadCubicSeed cubicPlane;
    std::cout << cubicPlane.infoString();

    
    std::cout << "TriHexSeed Memory Requirements: " << std::endl;
    std::cout << std::setw(20) << "refinement level" << std::setw(20) << "nverts" 
        << std::setw(20) << "nedges" << std::setw(20) << "nfaces" << std::endl;
    for (int i=0; i<maxRec; ++i) {
        triPlane.determineMaxAllocations(nv, ne, nf, i);
        std::cout << std::setw(20) << i << std::setw(20) << nv << std::setw(20) << ne  
            << std::setw(20) << nf << std::endl;
    }
    
    std::cout << "QuadRectSeed Memory Requirements: " << std::endl;
    std::cout << std::setw(20) << "refinement level" << std::setw(20) << "nverts" 
        << std::setw(20) << "nedges" << std::setw(20) << "nfaces" << std::endl;
    for (int i=0; i<maxRec; ++i) {
        quadPlane.determineMaxAllocations(nv, ne, nf, i);
        std::cout << std::setw(20) << i << std::setw(20) << nv << std::setw(20) << ne  
            << std::setw(20) << nf << std::endl;
    }

    std::cout << "QuadCubicSeed Memory Requirements:" << std::endl;
    std::cout << std::setw(20) << "refinement level" << std::setw(20) << "nverts" 
        << std::setw(20) << "nedges" << std::setw(20) << "nfaces" << std::endl;
    for (int i=0; i<maxRec; ++i) {
        cubicPlane.determineMaxAllocations(nv, ne, nf, i);
        std::cout << std::setw(20) << i << std::setw(20) << nv << std::setw(20) << ne  
            << std::setw(20) << nf << std::endl;
    }

    std::cout << "IcosTriSphereSeed Memory Requirements: " << std::endl;
    std::cout << std::setw(20) << "refinement level" << std::setw(20) << "nverts" 
        << std::setw(20) << "nedges" << std::setw(20) << "nfaces" << std::endl;
    for (int i=0; i<maxRec; ++i) {
        icosTriSphere.determineMaxAllocations(nv, ne, nf, i);
        std::cout << std::setw(20) << i << std::setw(20) << nv << std::setw(20) << ne  
            << std::setw(20) << nf << std::endl;
    }
    
    std::cout << "CubedSphereSeed Memory Requirements: " << std::endl;
    std::cout << std::setw(20) << "refinement level" << std::setw(20) << "nverts" 
        << std::setw(20) << "nedges" << std::setw(20) << "nfaces" << std::endl;
    for (int i=0; i<maxRec; ++i) {
        cubedSphere.determineMaxAllocations(nv, ne, nf, i);
        std::cout << std::setw(20) << i << std::setw(20) << nv << std::setw(20) << ne  
            << std::setw(20) << nf << std::endl;
    }


return 0;
}

