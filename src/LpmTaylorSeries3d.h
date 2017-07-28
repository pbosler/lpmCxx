#ifndef _LPM_TAYLOR_SERIES_3D_H_
#define _LPM_TAYLOR_SERIES_3D_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmMultiIndex.h"
#include "LpmCoords.h"
#include "LpmField.h"
#include "LpmFaces.h"
#include <map>
#include <vector>
#include <string>

namespace Lpm {

class TaylorSeries3d {
    public:
        virtual ~TaylorSeries3d() {};

        inline scalar_type getCoeff(const MultiIndex& k) const {return coeffs.at(k);}
        
        inline scalar_type getMoment(const MultiIndex& k) const {return moments.at(k);}
        
        scalar_type sum() const;
        
        std::string infoString() const;
        
        virtual void computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param = 0.0) = 0;
        
        virtual std::unique_ptr<TaylorSeries3d> createEmpty() const = 0;
        
        void computeMoments(const std::shared_ptr<Coords> crds, const std::vector<index_type>& crdInds, 
            const XyzVector& cntd, const std::shared_ptr<Field> srcStr);
            
        void computeMoments(const std::shared_ptr<Coords> crds, const std::vector<index_type>& crdInds, 
            const XyzVector& cntd, const std::shared_ptr<Field> srcVals, const std::shared_ptr<Field> srcWeights);
        
        void computeMoments(const std::shared_ptr<Faces> faces, const std::vector<index_type>& crdInds, 
            const XyzVector& cntd, const std::shared_ptr<Field> srcVals);
        
    protected:    
        TaylorSeries3d(const int maxSeriesOrder);
        std::map<MultiIndex, scalar_type> coeffs;
        std::map<MultiIndex, scalar_type> moments;
        int maxP;
        
        
};

class SphereGreensSeries : public TaylorSeries3d {
    public:
        SphereGreensSeries(const int maxSeriesOrder) : TaylorSeries3d(maxSeriesOrder) {};
        
        inline std::unique_ptr<TaylorSeries3d> createEmpty() const {
            return std::unique_ptr<TaylorSeries3d>(new SphereGreensSeries(maxP));
        }
        
        void computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param = 0.0);
};

class SecondOrderDelta3dSeries : public TaylorSeries3d {
    public:
        SecondOrderDelta3dSeries() : TaylorSeries3d(2) {};
        
        inline std::unique_ptr<TaylorSeries3d> createEmpty() const 
            {return std::unique_ptr<TaylorSeries3d>( new SecondOrderDelta3dSeries());}
        
        void computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param = 0.01);
};

class SphereDeltaSeries : public TaylorSeries3d {
    public:
        SphereDeltaSeries() : TaylorSeries3d(2) {};
        
        inline std::unique_ptr<TaylorSeries3d> createEmpty() const {
            return std::unique_ptr<TaylorSeries3d>(new SphereDeltaSeries());}
        
        void computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param = 0.01);
};

}

#endif
