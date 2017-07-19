#ifndef _LPM_TREE_SUM_H_
#define _LPM_TREE_SUM_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOctree.h"
#include "LpmXyzVector.h"
#include "LpmCoords.h"
#include "LpmField.h"
#include <vector>
#include <cmath>

namespace Lpm {

struct MultiIndex {
    MultiIndex(const unsigned i1 = 0, const unsigned i2 = 0, const unsigned i3 = 0) : k1(i1), k2(i2), k3(i3) {};

    unsigned k1;
    unsigned k2;
    unsigned k3;
    
    inline unsigned magnitude() const {return k1 + k2 + k3;}
    
    XyzVector vectorPower(const XyzVector& vec) const;
};

bool operator < (const MultiIndex& left, const MultiIndex& right);
bool operator > (const MultiIndex& left, const MultiIndex& right);
bool operator >= (const MultiIndex& left, const MultiIndex& right);
bool operator <= (const MultiIndex& left, const MultiIndex& right);
bool operator == (const MultiIndex& left, const MultiIndex& right);
bool operator != (const MultiIndex& left, const MultiIndex& right);


struct TaylorCoeffs {
        TaylorCoeffs(const int maxSeriesOrder);
        
//         inline scalar_type getCoeff(const int k1, const int k2, const int k3) const 
//             {return _coeffs.at(MultiIndex(k1,k2,k3));}
//         inline scalar_type getCoeff(const MultiIndex& kInd) const {return _coeffs.at(kInd);}
    
        void computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcCentroid) = 0;
        
        std::map<MultiIndex, scalar_type> vals;
        
//         inline index_type coeffCount(const int p) const {return (p*p*p + 6 * p*p + 11 * p + 6)/6;}
};


struct TreeSumNode : public TreeNode {
    TreeSumNode();
    TreeSumNode(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio, const int maxSeriesOrder, 
        const std::shared_ptr<Field> srcStrength);
    TreeSumNode(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio, const int maxSeriesOrder, 
        const std::shared_ptr<Field> srcVals, const std::shared_ptr<Field> srcWeights);
    
    void computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcStrength);
    void computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcVals, 
        const std::shared_ptr<Field> srcWeights);
    
    inline bool isFar(const XyzVector& tgtVec, const int maxTreeDepth, const scalar_type nuParam) const {
        const scalar_type hnu = std::pow(std::pow(2.0, -maxTreeDepth), nuParam);
        const XyzVector vec = tgtVec - box.centroid();
        return box.radius <= hnu * vec.magnitude();
    }
    
    inline bool isNear(const XyzVector& tgtVec, const int maxTreeDepth, const scalar_type nuParam) const {
        return !isFar(tgtVec, maxTreeDepth, nuParam);
    }

    std::map<MultiIndex, scalar_type> momentsA;
    std::map<MultiIndex, scalar_type> momentsB;
    std::map<MultiIndex, scalar_type> momentsC;
    
    TaylorCoeffs coeffs;
};

}

#endif
