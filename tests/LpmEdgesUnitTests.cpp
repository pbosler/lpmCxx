#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmEuclideanCoords.h"
#include "LpmEdges.h"
#include <memory>

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
	
	
	Edges edges(20, ec, lagec);
	edges.insert(0, 1, 0, -1); // 0
	edges.insert(1, 2, 1, -1); // 1
 	edges.insert(2, 3, 1, -1); // 2
	edges.insert(3, 4, 2, -1); // 3
	edges.insert(4, 5, 2, -1); // 4
	edges.insert(5, 6, 3, -1); // 5
	edges.insert(6, 7, 3, -1); // 6
	edges.insert(7, 0, 0, -1); // 7
	edges.insert(1, 8, 0, 1); // 8
	edges.insert(8, 5, 3, 2); // 9
	edges.insert(3, 8, 1, 2); // 10
	edges.insert(8, 7, 0, 3); // 11
	
	std::cout << "created " << edges.n() << " edges using space allocated for " << edges.nMax() << " edges." << std::endl;
	
	for (int i = 0; i < edges.n(); ++i)
	    std::cout << "edge " << i << " is on boundary : "  << (edges.onBoundary(i) ? "true" : "false") << std::endl;
    std::cout << "edge 0 record: " << edges.orig(0) << ", " << edges.dest(0) << ", " << 
        edges.leftFace(0) << ", " << edges.rightFace(0) << std::endl;
    std::cout << "\tedge 0 coords: origin " << edges.origCoord(0) << " destination " << edges.destCoord(0) << 
        " edgeVector " << edges.edgeVector(0) << std::endl;
    std::cout << "edge 0 midpoint should be (-1, 0.5). Computed midpoint is " << edges.midpoint(0) << std::endl;
    
    edges.divide(11);
    std::pair<index_type, index_type> kids = edges.children(11);
    std::cout << "divided edge 11.  Its children are: " << kids.first << " and " << kids.second << std::endl;
    std::cout << "\tadded new coordinate " << ec->getVec(ec->n()-1) << std::endl;
    lagcrdstr = lagec->listAllCoords();
	std::cout << lagcrdstr << std::endl;
    for (int i = 0; i < edges.n(); ++i)
	    std::cout << "edge " << i << " is on boundary : "  << (edges.onBoundary(i) ? "true" : "false") << std::endl;
return 0;
}
