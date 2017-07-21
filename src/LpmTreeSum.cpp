#include "LpmTreeSum.h"
#include "LpmMultiIndex.h"
#include "LpmTaylorSeries3d.h"
#include "LpmOutputMessage.h"
#include <cmath>
#include <typeinfo>
#include <limits>
#include <exception>

namespace Lpm {

TreeSumNode::TreeSumNode(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio, 
    const ScalarKernel& kernel, const int maxSeriesOrder) :
    Treenode(crds, maxAspectRatio) {
    
    const scalar_type fill_num = std::numeric_limits<scalar_type>::max();
    
    if (typeid(kernel) == typeid(SphereGreensFn) ) {
        coeffs = std::unique_ptr<TaylorCoeffs>(new SphereGreensCoeffs(maxSeriesOrder));
        
        scalarMoments.emplace(MultiIndex(0, 0, 0), fill_num);
        
        for (int k1 = 1; k1 <= maxSeriesOrder; ++k1) {
            scalarMoments.emplace(MultiIndex(k1, 0, 0), fill_num);
        }
        for (int k2 = 1; k2 <= maxSeriesOrder; ++k2) {
            scalarMoments.emplace(MultiIndex(0, k2, 0), fill_num);
        }
        for (int k1 = 1; k1 <= maxSeriesOrder; ++k1) {
            for(int k2 = 1; k2 <= maxSeriesOrder - k1; ++k2) {
                scalarMoments.emplace(MultiIndex(k1, k2, 0), fill_num);
            }
        }
        
        for (int k3 = 1; k3 <= maxSeriesOrder; ++k3) {
            scalarMoments.emplace(MultiIndex(0, 0, k3), fill_num);
            for (int k1 = 1; k1 <= maxSeriesOrder - k3; ++k1) {
                scalarMoments.emplace(MultiIndex(k1, 0, k3), fill_num);
            }
            for (int k2 = 1; k2 <= maxSeriesOrder - k3; ++k2) {
                scalarMoments.emplace(MultiIndex(0, k2, k3), fill_num);
            }
            for (int k1 = 1; k1 <= maxSeriesOrder - k3; ++k1) {
                for (int k2 = 1; k2 <= maxSeriesOrder - k3 - k1; ++k2) {
                    scalarMoments.emplace(MultiIndex(k1, k2, k3), fill_num);
                }
            }
        }
    }
    else {
        OutputMessage errMsg("Unrecognized kernel type.", OutputMessage::errorPriority, "TreeSumNode::TreeSumNode");
        log->logMessage(errMsg);
        throw std::invalid_argument("unrecognized kernel type.");
    }
};

void TreeSumNode::computeCoeffs(const XyzVector& tgtVec, const scalar_type param) {
    const XyzVector cntd = box.centroid();
    coeffs->computeCoeffs(tgtVec, cntd, param);
}

void TreeSumNode::computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcStrength) {
    const XyzVector cntd = box.centroid();
    const scalar_type moment = 0.0;
    for (index_type i = 0; i < nCoords(); ++i) {
        const XyzVector dVec = crds->getVec(coordsContained[i]) - cntd;
    }
}


}