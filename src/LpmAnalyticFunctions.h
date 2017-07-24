#ifndef _LPM_ANALYTIC_FUNCTION_H_
#define _LPM_ANALYTIC_FUNCTION_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include <cmath>

namespace Lpm {

class AnalyticFunction {
    public:
        inline virtual scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const 
            {return 0.0;}
        inline virtual scalar_type evaluateScalar(const XyzVector& crdVec) const {return 0.0;}
        inline virtual XyzVector evaluateVector(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const 
            {return XyzVector();}
        inline virtual XyzVector evaluateVector(const XyzVector& crdVec) const {return XyzVector();}
};

class SineWave3D : public AnalyticFunction {
    public:
        SineWave3D(const int kk, const int ll, const int mm) : k(kk), l(ll), m(mm) {};
    
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const override {
            return std::sin(2.0 * PI * k * x) * std::sin(2.0 * PI * l * y) * std::sin(2.0 * PI * m * z);
        };
        scalar_type evaluateScalar(const XyzVector& crdVec) const override {
            return evaluateScalar(crdVec.x, crdVec.y, crdVec.z);
        };
    protected:
        
        int k;
        int l;
        int m;
};

class Gaussian3D : public AnalyticFunction {
    public:
        Gaussian3D(const scalar_type b = 1.0) : _b(b) {};
    
        inline scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const 
            override {
            return std::exp(-_b * _b * ( x * x + y* y + z*z ));
        };
        inline scalar_type evaluateScalar(const XyzVector& crdVec) const override {
            return evaluateScalar(crdVec.x, crdVec.y, crdVec.z);
        };
        
    protected:        
        scalar_type _b;
};

class radial2dsource : public AnalyticFunction {
    public:
        radial2dsource() {};
    
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            const scalar_type rr = std::sqrt(x*x + y*y);
            const scalar_type rrdenom = rr / (ZERO_TOL * ZERO_TOL + rr * rr);
            return (std::abs(rr) <= 1.0 ? -0.5 * PI * (std::sin(PI * rr) * rrdenom + PI * std::cos(PI * rr)) : 0.0);
        };
        scalar_type evaluateScalar(const XyzVector& crdVec) const {
            return evaluateScalar(crdVec.x, crdVec.y);
        };
        XyzVector evaluateVector(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
        XyzVector evaluateVector(const XyzVector& crdVec) const {
            return XyzVector(0.0, 0.0, 0.0);
        }
}; 

class exact2dpotential : public AnalyticFunction {
    public:
        exact2dpotential() {};
        
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const override {
            const scalar_type rr = std::sqrt(x*x + y*y);
            return (std::abs(rr) <= 1.0 ? 0.5 * (1.0 + std::cos(PI * rr)) : 0.0);
        };
        scalar_type evaluateScalar(const XyzVector& crdVec) const override{
            return evaluateScalar(crdVec.x, crdVec.y);
        };
};

class sphereHarmonic54 : public AnalyticFunction {
    public:
        sphereHarmonic54() {}
        
        scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const override {
            const scalar_type lon = longitude(x, y);
            return 30.0 * std::cos(4.0 * lon) * legendre54(z);
        }
        
        scalar_type evaluateScalar(const XyzVector& crdVec) const override {
            return evaluateScalar(crdVec.x, crdVec.y, crdVec.z);
        }
        
        inline scalar_type legendre54(const scalar_type& z) const {
            return z * (z * z - 1.0) * (z * z - 1.0);
        }
};

}

#endif