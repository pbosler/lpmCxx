#ifndef _LPM_VECTOR_KERNEL_H_
#define _LPM_VECTOR_KERNEL_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include <cmath>

namespace Lpm {

class VectorKernel {
    public:
        virtual ~VectorKernel() {};
        virtual XyzVector evaluate(const XyzVector& tgtLoc, const XyzVector& srcLoc) const = 0;
        virtual bool isSingular() const = 0;
};

class BiotSavartSphere : public VectorKernel {
    public:
        BiotSavartSphere(const scalar_type delta = 0.0, const scalar_type sphRadius = 1.0) : _delta(delta),
            _sphRadius(sphRadius) {};
        
        inline XyzVector evaluate(const XyzVector& tgtLoc, const XyzVector& srcLoc) const {
            XyzVector result = tgtLoc.crossProduct(srcLoc);
            const scalar_type denom = 1.0 / (_sphRadius * _sphRadius - tgtLoc.dotProduct(srcLoc) + _delta * _delta);
            result.scale(denom);
            return result;
        }
        
        inline scalar_type sumMultiplier() const {
            return -1.0 / (4.0 * PI * _sphRadius);
        }
        
        inline bool isSingular() const {return _delta == 0.0; }
        
    protected:
        scalar_type _delta;
        scalar_type _sphRadius;
};

}

#endif
