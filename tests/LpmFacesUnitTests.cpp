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
    ec->insert(0.0, 0.0);
    ec->insert(0.5, 0.866025403784438597);
    ec->insert(-0.5, 0.866025403784438597);
    ec->insert(-1.0, 0.0);
    ec->insert(-0.5, -0.866025403784438597);
    ec->insert(0.5, -0.866025403784438597);
    ec->insert(1.0, 0.0);
    ec->insert(0.0, 0.577350269189625731);
    ec->insert(-0.5, 0.288675134594812810);
    ec->insert(-0.5, -0.288675134594812810);
    ec->insert(0.0, -0.577350269189625731);
    ec->insert(0.5, -0.288675134594812810);
    ec->insert(0.5, 0.288675134594812810);

	std::cout << "created " << ec->n() << " coordinates in the plan; space allocated for " << ec->nMax() << std::endl;
	
	
	std::shared_ptr<EuclideanCoords>  lagec(new EuclideanCoords(*ec));
	std::cout << "created " << lagec->n() << " lagrangian coordinates using copy constructor." << std::endl;
	std::string lagcrdstr = lagec->listAllCoords();
	std::cout << lagcrdstr << std::endl;
	
	std::shared_ptr<Edges> edges(new Edges(20, ec, lagec));
	edges->insert(0, 1, 0, 5); // 0
    edges->insert(1, 2, 0, -1); // 1
    edges->insert(2, 3, 1, -1); // 2
    edges->insert(3, 4, 2, -1); // 3 
    edges->insert(4, 5, 3, -1); // 4
    edges->insert(5, 6, 4, -1); // 5
    edges->insert(6, 1, 5, -1); // 6
    edges->insert(2, 0, 0, 1); // 7
    edges->insert(3, 0, 1, 2); // 8
    edges->insert(4, 0, 2, 3); // 9
    edges->insert(0, 5, 4, 3); // 10
    edges->insert(0, 6, 5, 4); // 11;
    
    std::cout << "created " << edges->n() << " edges\n";
    
// 0		1		0			5 // 0
// 1		2		0			-1 // 1
// 2		3		1			-1 // 2
// 3		4		2			-1 // 3
// 4		5		3			-1 // 4
// 5		6		4			-1 // 5
// 6		1		5			-1 // 6
// 2		0		0			1 // 7
// 3		0		1			2 // 8
// 4		0		2			3 // 9
// 0		5		4			3 // 10
// 0		6		5			4 // 11
	
	std::shared_ptr<Faces> faces(new TriFaces(20, edges, ec, lagec));
	std::vector<index_type> f0edges = {0, 1, 7};
	std::vector<index_type> f1edges = {2, 8, 7};
	std::vector<index_type> f2edges = {9, 8, 3};
	std::vector<index_type> f3edges = {9, 4, 10};
	std::vector<index_type> f4edges = {5, 11, 10};
	std::vector<index_type> f5edges = {0, 11, 6};
    faces->insert(f0edges);
    faces->insert(f1edges);
    faces->insert(f2edges);
    faces->insert(f3edges);
    faces->insert(f4edges);
    faces->insert(f5edges);
    
    std::cout << "created " << faces->n() << " faces\n";
    
    const scalar_type exactArea = 6 * (0.5 * 0.866025403784438597);
    std::cout << "exact surfArea = " << exactArea << std::endl;
    for (int i = 0; i < faces->n(); ++i)
        std::cout << i << ": " << faces->area(i) << std::endl;
    
    
    faces->resetAreas();
    for (int i = 0; i < faces->n(); ++i)
        std::cout << i << ": " << faces->area(i) << std::endl;
    
    faces->divide(0);
    
    std::cout << "faces->n() = " << faces->n() << "; faces->area(i) = " << std::endl;
    for (int i = 0; i < faces->n(); ++i)
        std::cout << i << ": " << faces->area(i) << std::endl;
    
    faces->resetAreas();
return 0;
}
