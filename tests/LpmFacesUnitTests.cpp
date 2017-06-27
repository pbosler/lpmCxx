#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmFaces.h"
#include "LpmEdges.h"
#include "LpmEuclideanCoords.h"
#include "LpmTriFaces.h"
//#include "LpmQuadFaces.h"
#include <memory>
#include <iostream>

using namespace Lpm;

int main(int argc, char* argv[]) {
    std::shared_ptr<EuclideanCoords> ec(new EuclideanCoords(20));
    ec->insert( -1.0, 1.0 ); //0
    ec->insert( -1.0, 0.0 ); //1
    ec->insert( -1.0,-1.0 ); //2
	ec->insert( 0.0, -1.0 ); //3
	ec->insert( 1.0, -1.0 ); //4: 1.0		-1.0 
	ec->insert(1.0, 0.0 ); // 5
	ec->insert(1.0, 1.0 ); //6
	ec->insert(0.0, 1.0 ); //7
	ec->insert(0.0, 0.0 ); //8
	ec->insert(-0.5, 0.5 ); //9
	XyzVector vec1(-0.5, -0.5 );
	XyzVector vec2(0.5,	-0.5 );
	XyzVector vec3(0.5,  0.5 );
	ec->insert(vec1); //10
	ec->insert(vec2); //11
	ec->insert(vec3); //12
	std::cout << "created " << ec->n() << " coordinates in the plan; space allocated for " << ec->nMax() << std::endl;
	
	
	std::shared_ptr<EuclideanCoords>  lagec(new EuclideanCoords(*ec));
	std::cout << "created " << lagec->n() << " lagrangian coordinates using copy constructor." << std::endl;
	std::string lagcrdstr = lagec->listAllCoords();
	std::cout << lagcrdstr << std::endl;
	
	std::shared_ptr<Edges> edges(new Edges(20, ec, lagec));
	edges->insert(0, 1, 0, -1); // 0
	edges->insert(1, 2, 1, -1); // 1
 	edges->insert(2, 3, 1, -1); // 2
	edges->insert(3, 4, 2, -1); // 3
	edges->insert(4, 5, 2, -1); // 4
	edges->insert(5, 6, 3, -1); // 5
	edges->insert(6, 7, 3, -1); // 6
	edges->insert(7, 0, 0, -1); // 7
	edges->insert(1, 8, 0, 1); // 8
	edges->insert(8, 5, 3, 2); // 9
	edges->insert(3, 8, 1, 2); // 10
	edges->insert(8, 7, 0, 3); // 11
	
	
return 0;
}
