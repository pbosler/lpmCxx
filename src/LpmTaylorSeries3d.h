#ifndef _LPM_TAYLOR_SERIES_3D_H_
#define _LPM_TAYLOR_SERIES_3D_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmMultiIndex.h"
#include <map>
#include <vector>
#include <string>

namespace Lpm {

class TaylorCoeffs {
    public:
        virtual ~TaylorCoeffs() {};

        inline scalar_type getCoeff(const MultiIndex& kk) const {return aa.at(kk);}
        
        std::string infoString() const;
        
        virtual void computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param = 0.0) = 0;
        
    protected:    
        TaylorCoeffs(const int maxSeriesOrder);
        std::map<MultiIndex, scalar_type> aa;
        int maxP;
        
        
};

class SphereGreensCoeffs : public TaylorCoeffs {
    public:
        SphereGreensCoeffs(const int maxSeriesOrder) : TaylorCoeffs(maxSeriesOrder) {};
        
        void computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param = 0.0);
};

}

#endif
