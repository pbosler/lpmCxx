#ifndef _LPM_TREE_SUM_H_
#define _LPM_TREE_SUM_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmLogger.h"
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

struct SumNode : public Node {
    typedef std::unique_ptr<TaylorSeries3d> series_ptr_type;

    SumNode(const Box3d& bbox, Node* pparent = NULL, const std::vector<index_type>& crdInds = std::vector<index_type>(),
         const int maxSeriesOrder = 0, ScalarKernel* kernel = NULL);

    series_ptr_type series;
    
    void computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcStrength);
    void computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcVals, 
        const std::shared_ptr<Field> srcWeights);
    void computeCoeffs(const XyzVector& tgtVec, const scalar_type param = 0.0);
    
    inline bool isFar(const XyzVector& tgtVec, const int maxTreeDepth, const scalar_type nuParam = 1.0) const {
        const scalar_type hnu = std::pow(std::pow(2.0, -maxTreeDepth), nuParam);
        const XyzVector vec = tgtVec - box.centroid();
        return box.radius() <= hnu * vec.magnitude();
    }
    
    inline bool isNear(const XyzVector& tgtVec, const int maxTreeDepth, const scalar_type nuParam) const {
        return !isFar(tgtVec, maxTreeDepth, nuParam);
    }
    
};

class TreeSum : public Tree {
    public:
        TreeSum(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio, 
            const std::shared_ptr<ScalarKernel> kernel, const int maxSeriesOrder, const scalar_type seriesParam = 0.0,
            const int prank = 0);

        void buildTree(const index_type maxCoordsPerNode) override;

    protected:
        std::weak_ptr<ScalarKernel> _kernel;
        int _maxP;
        scalar_type _seriesParam;
    
        void generateTree(Node* node, const index_type maxCoordsPerNode) override;
};


}

#endif
