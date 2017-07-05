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
    
    protected:
};

class PlanarGreensFnFreeBoundaries : public ScalarKernel {
    public:
        PlanarGreensFnFreeBoundaries() {};
        
        inline scalar_type evaluate(const XyzVector& tgtLoc, const XyzVector &srcLoc) const {
            const XyzVector interactionVector = tgtLoc - srcLoc;
            return std::log( interactionVector.magnitudeSquared()) / (4.0 * PI);
        }
};

}

#endif
