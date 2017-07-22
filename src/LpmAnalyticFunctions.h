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

}

#endif