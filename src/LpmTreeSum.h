#ifndef _LPM_TREE_SUM_H_
#define _LPM_TREE_SUM_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOctree.h"
#include "LpmXyzVector.h"
#include "LpmCoords.h"
#include "LpmField.h"
#include "LpmMultiIndex.h"
#include "LpmTaylorSeries3d.h"
#include "LpmScalarKernel.h"
#include <vector>
#include <cmath>
#include <memory>

namespace Lpm {

struct TreeSumNode : public Treenode {
    TreeSumNode(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio, 
        const ScalarKernel& kernel, const int maxSeriesOrder, const scalar_type seriesParam = 0.0);
        
    TreeSumNode(const Box3d& bbox, const std::shared_ptr<TreeSumNode> pparent = NULL, 
        const std::vector<index_type>& crdInds = std::vector<index_type>(), const scalar_type maxAspectRatio = 1.0);
        
    void computeCoeffs(const XyzVector& tgtVec, const scalar_type param = 0.0);
    
    void computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcStrength);
    void computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcVals, 
        const std::shared_ptr<Field> srcWeights);
    
    inline bool isFar(const XyzVector& tgtVec, const int maxTreeDepth, const scalar_type nuParam = 1.0) const {
        const scalar_type hnu = std::pow(std::pow(2.0, -maxTreeDepth), nuParam);
        const XyzVector vec = tgtVec - box.centroid();
        return box.radius() <= hnu * vec.magnitude();
    }
    
    inline bool isNear(const XyzVector& tgtVec, const int maxTreeDepth, const scalar_type nuParam) const {
        return !isFar(tgtVec, maxTreeDepth, nuParam);
    }

    std::unique_ptr<TaylorSeries3d> series;
};

void generateTree(std::shared_ptr<TreeSumNode> node, std::shared_ptr<Coords> crds, const index_type maxCoordsPerNode);

void computeAllMoments(std::shared_ptr<TreeSumNode> tree, const std::shared_ptr<Coords> crds, 
    const std::shared_ptr<Field> srcStrength);

void computeAllMoments(std::shared_ptr<TreeSumNode> tree, const std::shared_ptr<Coords> crds, 
    const std::shared_ptr<Field> srcVals, const std::shared_ptr<Field> srcWeights);

}

#endif
