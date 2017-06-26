#ifndef _LPM_ANALYTIC_FUNCTION_H_
#define _LPM_ANALYTIC_FUNCTION_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"

namespace Lpm {

class AnalyticFunction {
    public:
        virtual scalar_type evaluateScalar(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const = 0;
        virtual scalar_type evaluateScalar(const XyzVector& crdVec) const = 0;
        virtual XyzVector evaluateVector(const scalar_type x, const scalar_type y, const scalar_type z = 0.0) const = 0;
        virtual XyzVector evaluateVector(const XyzVector& crdVec) const = 0;
};


}

#endif