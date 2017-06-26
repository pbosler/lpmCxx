#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmField.h"
#include "LpmCoords.h"
#include "LpmEuclideanCoords.h"
#include "LpmAnalyticFunctions.h"
#include "LpmXyzVector.h"
#include <cmath>

using namespace Lpm;

class SineWave3D : public AnalyticFunction {
    public:
        SineWave3D(const int kk, const int ll, const int mm) : k(kk), l(ll), m(mm) {};
    
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            return std::sin(2.0 * PI * k * x) * std::sin(2.0 * PI * l * y) * std::sin(2.0 * PI * m * z);
        };
        scalar_type evaluateScalar(const XyzVector& crdVec) const {
            return evaluateScalar(crdVec.x, crdVec.y, crdVec.z);
        };
        XyzVector evaluateVector(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
        XyzVector evaluateVector(const XyzVector& crdVec) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
    protected:
        
        int k;
        int l;
        int m;
};

class Gaussian3D : public AnalyticFunction {
    public:
        Gaussian3D(const scalar_type b = 1.0) : _b(b) {};
    
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            return std::exp(-_b * _b * ( x * x + y* y + z*z ));;
        };
        scalar_type evaluateScalar(const XyzVector& crdVec) const {
            return evaluateScalar(crdVec.x, crdVec.y, crdVec.z);
        };
        XyzVector evaluateVector(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
        XyzVector evaluateVector(const XyzVector& crdVec) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
    protected:
        
        scalar_type _b;
};

int main (int argc, char* argv[] ) {

    const int n = 13;

    Field sineScalar(n, 1, "f(x)");
    Field gaussScalar(n, 1, "g(x)");
    const SineWave3D sw3(2, 5, 7);
    const Gaussian3D gs3(1.0);
    
    EuclideanCoords ec(13);
    
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


