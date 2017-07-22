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
    const ScalarKernel& kernel, const int maxSeriesOrder, const scalar_type seriesParam) : 
    Treenode(crds, maxAspectRatio) {
    
    if (typeid(kernel) == typeid(SphereGreensFn)) {
        series = std::unique_ptr<TaylorSeries3d>(new SphereGreensSeries(maxSeriesOrder));
    }
    else if (typeid(kernel) == typeid(SecondOrderDelta3d)) {
        series = std::unique_ptr<TaylorSeries3d>(new SecondOrderDelta3dSeries());
    }
    else {
        OutputMessage errMsg("Unrecognized kernel", OutputMessage::errorPriority, "TreeSumNode::TreeSumNode");
        log->logMessage(errMsg);
        throw std::invalid_argument("Unrecognized kernel");
    }
};

TreeSumNode::TreeSumNode(const Box3d& bbox, const std::shared_ptr<TreeSumNode> pparent, 
    const std::vector<index_type>& crdInds, const scalar_type maxAspectRatio ) : Treenode(bbox, pparent, crdInds, maxAspectRatio) {
    if (pparent.use_count() > 0) {
        series = pparent->series->createEmpty();
    }
}

void TreeSumNode::computeCoeffs(const XyzVector& tgtVec, const scalar_type param) {
    const XyzVector cntd = box.centroid();
    series->computeCoeffs(tgtVec, cntd, param);
}



}