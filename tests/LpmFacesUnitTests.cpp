#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmFaces.h"
#include "LpmEdges.h"
#include "LpmEuclideanCoords.h"
#include "LpmTriFaces.h"
#include "LpmQuadFaces.h"
#include <memory>
#include <iostream>

using namespace Lpm;

int main(int argc, char* argv[]) {

    { // Triangular faces
        std::cout << "************** Test 1: Triangular Faces:\n";
        
        std::shared_ptr<EuclideanCoords> ec(new EuclideanCoords(20));
        ec->insert(0.0, 0.0); //0
        ec->insert(0.5, 0.866025403784438597); //1
        ec->insert(-0.5, 0.866025403784438597); //2
        ec->insert(-1.0, 0.0); //3
        ec->insert(-0.5, -0.866025403784438597);//4
        ec->insert(0.5, -0.866025403784438597);//5
        ec->insert(1.0, 0.0);//6
        ec->insert(0.0, 0.577350269189625731);//7
        ec->insert(-0.5, 0.288675134594812810);
        ec->insert(-0.5, -0.288675134594812810);
        ec->insert(0.0, -0.577350269189625731);
        ec->insert(0.5, -0.288675134594812810);
        ec->insert(0.5, 0.288675134594812810);

        std::cout << "created " << ec->n() << " coordinates in the plane; space allocated for " << ec->nMax() << std::endl;
    
    
        std::shared_ptr<EuclideanCoords>  lagec(new EuclideanCoords(*ec));
        std::cout << "created " << lagec->n() << " lagrangian coordinates using copy constructor." << std::endl;
        std::string lagcrdstr = lagec->listAllCoords();
        std::cout << lagcrdstr << std::endl;
    
        std::shared_ptr<Edges> edges(new Edges(50, ec, lagec));
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
        std::cout << "coords.n(), coords.nMax() = " << ec->n() << ", " << ec->nMax() << std::endl;
        std::cout << "edges.n(), edges.nMax() = " << edges->n() << ", " << edges->nMax() << std::endl;
        std::cout << "faces.n(), faces.nMax() = " << faces->n() << ", " << faces->nMax() << std::endl;
    
        std::cout << "faces->n() = " << faces->n() << "; faces->area(i) = " << std::endl;
        for (int i = 0; i < faces->n(); ++i)
            std::cout << i << ": " << faces->area(i) << std::endl;
    
        faces->resetAreas();
        
        faces->divide(8);
        std::cout << "coords.n(), coords.nMax() = " << ec->n() << ", " << ec->nMax() << std::endl;
        std::cout << "edges.n(), edges.nMax() = " << edges->n() << ", " << edges->nMax() << std::endl;
        std::cout << "faces.n(), faces.nMax() = " << faces->n() << ", " << faces->nMax() << std::endl;
        std::cout << "faces->n() = " << faces->n() << "; faces->area(i) = " << std::endl;
        for (int i = 0; i < faces->n(); ++i)
            std::cout << i << ": " << faces->area(i) << std::endl;
        faces->resetAreas();
    }
    {// Quadrilateral faces
        std::cout << "************** Test 2: Quadrilateral Faces:\n";
    
        std::shared_ptr<EuclideanCoords> ec(new EuclideanCoords(30));
        
        ec->insert(-1.0, 1.0 );
        ec->insert(-1.0, 0.0 );
        ec->insert(-1.0, -1.0 );
        ec->insert(0.0,	-1.0 );
        ec->insert(1.0,	-1.0 );
        ec->insert(1.0,	 0.0 );
        ec->insert(1.0,	 1.0 );
        ec->insert(0.0,	 1.0 );
        ec->insert(0.0,	 0.0 );
        ec->insert(-0.5, 0.5 );
        ec->insert(-0.5, -0.5 );
        ec->insert(0.5,	-0.5 );
        ec->insert(0.5,	 0.5 );
        
        std::cout << "created " << ec->n() << " coordinates in the plane; space allocated for " << ec->nMax() << std::endl;
        
        std::shared_ptr<Edges> edges(new Edges(50, ec));
        edges->insert(0, 1,	0, -1);
        edges->insert(1, 2,	1, -1);
        edges->insert(2, 3,	1, -1);
        edges->insert(3, 4,	2, -1);
        edges->insert(4, 5,	2, -1);
        edges->insert(5, 6,	3, -1);
        edges->insert(6, 7,	3, -1);
        edges->insert(7, 0,	0, -1);
        edges->insert(1, 8,	0, 1);
        edges->insert(8, 5,	3, 2);
        edges->insert(3, 8,	1, 2);
        edges->insert(8, 7, 0, 3);
        
        std::cout << "created " << edges->n() << " edges\n";
        
        std::shared_ptr<Faces> faces(new QuadFaces(20, edges, ec));
        
        std::vector<index_type> f0edges = {0,	8,	11, 7};
        std::vector<index_type> f1edges = {1,	2,	10,	8};
        std::vector<index_type> f2edges = {10,	3,	4,	9};
        std::vector<index_type> f3edges = {11,	9,	5,	6};
        
        faces->insert(f0edges);
        faces->insert(f1edges);
        faces->insert(f2edges);
        faces->insert(f3edges);
        
        std::cout << "created " << faces->n() << " faces\n";
        
    
        const scalar_type exactArea = 4.0;
        std::cout << "exact surfArea = " << exactArea << std::endl;
        for (int i = 0; i < faces->n(); ++i)
            std::cout << i << ": " << faces->area(i) << std::endl;
    
    
        faces->resetAreas();
        for (int i = 0; i < faces->n(); ++i)
            std::cout << i << ": " << faces->area(i) << std::endl;
    
        faces->divide(0);
        std::cout << "coords.n(), coords.nMax() = " << ec->n() << ", " << ec->nMax() << std::endl;
        std::cout << "edges.n(), edges.nMax() = " << edges->n() << ", " << edges->nMax() << std::endl;
        std::cout << "faces.n(), faces.nMax() = " << faces->n() << ", " << faces->nMax() << std::endl;
        std::cout << "faces->n() = " << faces->n() << "; faces->area(i) = " << std::endl;
        for (int i = 0; i < faces->n(); ++i)
            std::cout << i << ": " << faces->area(i) << std::endl;
    
        std::cout << "faces->n() = " << faces->n() << "; faces->area(i) = " << std::endl;
        for (int i = 0; i < faces->n(); ++i)
            std::cout << i << ": " << faces->area(i) << std::endl;
    
        faces->resetAreas();
        
        faces->divide(6);
        std::cout << "coords.n(), coords.nMax() = " << ec->n() << ", " << ec->nMax() << std::endl;
        std::cout << "edges.n(), edges.nMax() = " << edges->n() << ", " << edges->nMax() << std::endl;
        std::cout << "faces.n(), faces.nMax() = " << faces->n() << ", " << faces->nMax() << std::endl;
        std::cout << "faces->n() = " << faces->n() << "; faces->area(i) = " << std::endl;
        for (int i = 0; i < faces->n(); ++i)
            std::cout << i << ": " << faces->area(i) << std::endl;
        faces->resetAreas();
    }
return 0;
}
