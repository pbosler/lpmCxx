#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmField.h"
#include "LpmCoords.h"
#include "LpmEuclideanCoords.h"
#include "LpmAnalyticFunctions.h"
#include "LpmXyzVector.h"
#include <cmath>

using namespace Lpm;



int main (int argc, char* argv[] ) {

    const int n = 13;

    Field sineScalar(n, 1, "f(x)");
    Field gaussScalar(n, 1, "g(x)");
    const SineWave3D sw3(2, 5, 7);
    const Gaussian3D gs3(1.0);
    
    EuclideanCoords ec(13, CARTESIAN_3D_GEOMETRY);
    
    ec.insert( -1.0, 1.0 ); //0
    ec.insert( -1.0, 0.0 ); //1
    ec.insert( -1.0,-1.0 ); //2
	ec.insert( 0.0, -1.0 ); //3
	ec.insert( 1.0, -1.0 ); //4: 1.0		-1.0 
	ec.insert(1.0, 0.0 ); // 5
	ec.insert(1.0, 1.0 ); //6
	ec.insert(0.0, 1.0 ); //7
	ec.insert(0.0, 0.0 ); //8
	ec.insert(-0.5, 0.5 ); //9
	XyzVector vec1(-0.5, -0.5 );
	XyzVector vec2(0.5,	-0.5 );
	XyzVector vec3(0.5,  0.5 );
	ec.insert(vec1); //10
	ec.insert(vec2); //11
	ec.insert(vec3); //12
    
    sineScalar.initializeToScalarFunction(&ec, &sw3);
    gaussScalar.initializeToScalarFunction(&ec, &gs3);
    
    std::cout << "sine wave values : x, f(x)" << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << i << ": " << ec.getVec(i) << " f(x) = " << sineScalar.getScalar(i) << std::endl;
    }
    std::cout << "Gaussian values : x, g(x) " << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << i << ": " << ec.getVec(i) << " g(x) = " << gaussScalar.getScalar(i) << std::endl;
    }
return 0;
}


