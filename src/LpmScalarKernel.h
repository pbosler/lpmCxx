#ifndef _LPM_HEADER_TEMPLATE_H_
#define _LPM_HEADER_TEMPLATE_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include <cmath>

namespace Lpm {

class ScalarKernel {
    public:
        virtual ~ScalarKernel() {};
        virtual scalar_type evaluate(const XyzVector& tgtLoc, const XyzVector& srcLoc) const = 0;
        virtual bool isSingular() const = 0;
    
};

class PlanarGreensFnFreeBoundaries : public ScalarKernel {
    public:
        PlanarGreensFnFreeBoundaries() {};
        
        inline scalar_type evaluate(const XyzVector& tgtLoc, const XyzVector &srcLoc) const {
            const XyzVector interactionVector = tgtLoc - srcLoc;
            return std::log( interactionVector.magnitudeSquared()) / (4.0 * PI);
        }
        
        inline bool isSingular() const {return true;}
};

class SphereGreensFn : public ScalarKernel {
    public:
        SphereGreensFn(const scalar_type r) : _radius(r) {};
    
        inline scalar_type evaluate(const XyzVector& tgtLoc, const XyzVector& srcLoc) const {
            const scalar_type inarg = _radius * _radius - srcLoc.dotProduct(tgtLoc);
            return -std::log(inarg) / (4.0 * PI * _radius);
        }
        
        inline bool isSingular() const {return true;}
        
    protected:
        scalar_type _radius;
    
};

class SecondOrderDelta3d : public ScalarKernel {
    public:
        SecondOrderDelta3d(const scalar_type eps) : _eps(eps) {};
        
        inline scalar_type evaluate(const XyzVector& tgtVec, const XyzVector& srcVec) const {
            const XyzVector dvec = tgtVec - srcVec;
            const scalar_type r = dvec.magnitude();
            const scalar_type expfactor = std::exp(-std::pow(r/_eps, 3));
            const scalar_type eps3 = std::pow(_eps, 3);
            return 3.0 * expfactor / (4.0 * PI * eps3);
        }
        
        inline bool isSingular() const {return false;}
        
    protected:
        scalar_type _eps;
};

class SphereDelta : public ScalarKernel {
    public:
        SphereDelta(const scalar_type eps) : _eps(eps) {};
    
        inline scalar_type evaluate(const XyzVector& tgtVec, const XyzVector& srcVec) const {
            static const scalar_type bkgrnd = 1.0 / (4.0 * PI);
            const scalar_type dotProd = tgtVec.dotProduct(srcVec);
            const scalar_type numer = 2.0 * _eps * _eps * dotProd - (1.0 - dotProd) * (1.0 - dotProd);
            const scalar_type denom = 4.0 * PI * (1.0 - dotProd + _eps * _eps) * (1.0 - dotProd + _eps * _eps);
            return numer / denom - bkgrnd;
        }    
        
        inline bool isSingular() const {return false;}
    protected:
        scalar_type _eps;
};

}

#endif
