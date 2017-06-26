#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using namespace Lpm;

int main() {

    XyzVector vec3( 0.0, 0.0, 4.0);
    cout << " threeD vector : \n" << vec3 ;
    cout << " has magnitude = " << vec3.magnitude() << endl;
    XyzVector vec32( 1.0, 0.0, 0.0);
    XyzVector cp = vec3.crossProduct(vec32);
    cout << " cross product vector : \n" << cp ;
    cout << " has magnitude = " << cp.magnitude() << endl;
    cout << "\n";
      
    cout << "testing copy : " << endl;
    XyzVector vec3copy(vec3);
    cout << vec3copy << endl;
    
    cout <<  "testing assignment : " << endl;
    vec3 = vec32;
    cout << vec3 << endl;
    
    XyzVector a(0.0, 0.0);
    XyzVector b(1.0, 0.0);
    XyzVector c(1.0, 1.0);
    XyzVector d(0.0, 1.0);
    
    std::vector< XyzVector > tri(3);
    tri[0] = a;
    tri[1] = b;
    tri[2] = c;
    std::vector< XyzVector > quad(4);
    quad[0] = a;
    quad[1] = b;
    quad[2] = c;
    quad[3] = d;
    
    cout << "triangle area should be 0.5; computed area is " << triArea(a, b, c ) << " .\n";
    
    cout << "quad centroid should be (0.5, 0.5); computed centroid is " << centroid(quad) << ".\n";
    cout << "test += operator is part of centroid function.\n";
    
    cout << "testing + operator : " << endl;
    XyzVector oneone = b + d;
    cout << "sum should be (1,1); sum is " << oneone << ".\n";
    
    cout << "testing - operator : " << endl;
    XyzVector orig = oneone - b - d;
    cout << "difference should be (0,0); difference is " << orig << ".\n";
    
    cout << "test == operator : " << endl;
    bool ans = (orig == a);
    cout << "answer should be true; answer is " << ans << endl;

return 0;
}