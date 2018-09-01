#include "LpmAosMeshSeed.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace Lpm::Aos;

int main(int argc, char* argv[]) {

    TriHexSeed triPlane;
    
    std::cout << "id string = " << triPlane.idString() << std::endl;
    
    std::cout << "initializing from file." << std::endl;
    triPlane.initFromFile();
    std::cout << "returned." << std::endl;
    std::cout << triPlane.infoString();
    
    QuadRectSeed quadPlane;
    quadPlane.initFromFile();
    std::cout << quadPlane.infoString();
    
    IcosTriSphereSeed icosTriSphere;
    icosTriSphere.initFromFile();
    std::cout << icosTriSphere.infoString();
    
    CubedSphereSeed cubedSphere;
    cubedSphere.initFromFile();
    std::cout << cubedSphere.infoString();

return 0;
}

