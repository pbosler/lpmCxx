#include "LpmTreeSum.h"
#include "LpmMultiIndex.h"
#include "LpmTaylorSeries3d.h"
#include "LpmOutputMessage.h"
#include <cmath>
#include <limits>
#include <exception>

namespace Lpm {

SumNode::SumNode(const Box3d& bbox, Node* pparent, const std::vector<index_type>& crdInds, const int maxSeriesOrder,
    ScalarKernel* kernel) : Node(bbox, pparent, crdInds), momentsReady(false) {
    
    PlanarGreensFnFreeBoundaries* plane_ptr = dynamic_cast<PlanarGreensFnFreeBoundaries*>(kernel);
    SphereGreensFn* sphere_green_ptr = dynamic_cast<SphereGreensFn*>(kernel);
    SecondOrderDelta3d* delta_3d_ptr = dynamic_cast<SecondOrderDelta3d*>(kernel);
    SphereDelta* sphere_delta_ptr = dynamic_cast<SphereDelta*>(kernel);
    
    if (plane_ptr) {
        OutputMessage errMsg("Not implemented yet.", OutputMessage::errorPriority, "SumNode::SumNode");
        log->logMessage(errMsg);
        throw std::runtime_error("Planar Greens function Taylor series not implemented yet.");
    }
    else if (sphere_green_ptr) {
        series = series_ptr_type(new SphereGreensSeries(maxSeriesOrder));
    }
    else if (delta_3d_ptr) {
        series = series_ptr_type(new SecondOrderDelta3dSeries());
    }
    else if (sphere_delta_ptr) {
        series = series_ptr_type(new SphereDeltaSeries());
    }
    else {
        OutputMessage errMsg("Cannot determine series type from kernel.", OutputMessage::errorPriority, "SumNode::SumNode");
        log->logMessage(errMsg);
        throw std::invalid_argument("Cannot determine series type from kernel.");
    }
}


void SumNode::computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcStrength) {
    series->computeMoments(crds, coordsContained, box.centroid(), srcStrength);
    momentsReady = true;
}

void SumNode::computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcVals, 
    const std::shared_ptr<Field> srcWeights) {
    series->computeMoments(crds, coordsContained, box.centroid(), srcVals, srcWeights);
    momentsReady = true;
}

void SumNode::computeCoeffs(const XyzVector& tgtVec, const scalar_type param){
    series->computeCoeffs(tgtVec, box.centroid(), param);
}
    
TreeSum::TreeSum(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio, 
    const std::shared_ptr<ScalarKernel> kernel, const int maxSeriesOrder, const scalar_type seriesParam,
    const int prank, const scalar_type farFieldParam) : Tree(crds, maxAspectRatio, prank), _maxP(maxSeriesOrder), _kernel(kernel), 
    _seriesParam(seriesParam), _nuParam(farFieldParam)
{
    std::unique_ptr<Node> new_root(new SumNode(_root->box, NULL, _root->coordsContained, _maxP, _kernel.lock().get()));
    _root = std::move(new_root);
}

void TreeSum::buildTree(const index_type maxCoordsPerNode) {
    generateTree(_root.get(), maxCoordsPerNode);
    _depth = computeTreeDepth(_root.get());
}

void TreeSum::generateTree(Node* node, const index_type maxCoordsPerNode) {
    if (node->coordsContained.size() <= maxCoordsPerNode) {
        return;
    }
    else {
        //
        //  determine box dimensions to split
        //
        std::shared_ptr<Coords> crd_ptr = _crds.lock();
        bool splitDims[3];
        int splitCount = 0;
        const scalar_type edgeThreshold = node->box.longestEdge() / _maxAspectRatio;
        for (int i = 0; i < 3; ++i) {
            if (node->box.edgeLength(i) >= edgeThreshold) {
                splitDims[i] = true;
                splitCount += 1;
            }
            else {
                splitDims[i] = false;
            }
        }
        //
        //  make child boxes
        //
        std::vector<Box3d> kidboxes = node->box.bisectAlongDims(splitDims);
        for (int i = 0; i < kidboxes.size(); ++i) {
            //
            //  find coordinates contained by child box
            //
            std::vector<index_type> kidcoords;
            kidcoords.reserve(node->coordsContained.size());
            for (index_type j = 0; j < node->coordsContained.size(); ++j) {
                if (kidboxes[i].containsPoint(crd_ptr->getVec(node->coordsContained[j]))) {
                    kidcoords.push_back(node->coordsContained[j]);
                }
            }
            if (!kidcoords.empty()) {
                kidcoords.shrink_to_fit();
                node->kids.push_back(std::unique_ptr<Node>(new SumNode(kidboxes[i], node, kidcoords, _maxP, _kernel.lock().get())));
            }
        }
        if (node->kids.empty()) {
            OutputMessage errMsg("All kids are empty, this shouldn't happen.", OutputMessage::errorPriority, "Treee::buildTree");
            log->logMessage(errMsg);
            return;
        }
        else {
            _nnodes += node->kids.size();
            for (int i = 0; i < node->kids.size(); ++i) {
                shrinkBox(node->kids[i].get());  
                generateTree(node->kids[i].get(), maxCoordsPerNode);
            }
        }
    }
}

void TreeSum::setRecomputeMomentsTrue(Node* node) {
    SumNode* sum_ptr = dynamic_cast<SumNode*>(node);
    if (sum_ptr) {
        sum_ptr->momentsReady = false;
        if (sum_ptr->hasKids()) {
            for (int i = 0; i < sum_ptr->kids.size(); ++i) {
                setRecomputeMomentsTrue(sum_ptr->kids[i].get());
            }
        }
    }
    else {
        OutputMessage errMsg("failed to cast Node* to SumNode*", OutputMessage::errorPriority, "Tree::setRecomputeMomentsTrue");
        log->logMessage(errMsg);
        throw std::logic_error("");
    }
}

scalar_type TreeSum::computeSum(const XyzVector& tgtLoc, const std::shared_ptr<Field> srcStrength) {

}

scalar_type TreeSum::computeSum(const XyzVector& tgtLoc, const std::shared_ptr<Field> srcVals, 
    const std::shared_ptr<Field> srcWeights) {
    SumNode* root_ptr = dynamic_cast<SumNode*>(_root.get());
    scalar_type result = 0.0;
    recursiveSum(result, tgtLoc, srcVals, srcWeights, root_ptr, _depth);
    return result;
}

void TreeSum::recursiveSum(scalar_type& sum, const XyzVector& tgtLoc, const std::shared_ptr<Field> srcVals, 
    const std::shared_ptr<Field> srcWeights, SumNode* node, const index_type& tree_depth) {
    if (node->isFar(tgtLoc, tree_depth, _nuParam)) {
        //
        //  add Taylor series approximation
        //
        if (!node->momentsReady) {
            node->computeMoments(_crds.lock(), srcVals, srcWeights);
        }
        node->computeCoeffs(tgtLoc, _seriesParam);
        sum += node->seriesSum();
    }
    else {
        if (node->isLeaf()) {
            //
            //  add direct sum 
            //
            std::shared_ptr<ScalarKernel> kernel_ptr = _kernel.lock();
            if (kernel_ptr->isSingular()) {
                for (index_type i = 0; i < node->coordsContained.size(); ++i) {
                    const XyzVector srcVec = _crds.lock()->getVec(node->coordsContained[i]);
                    if (_crds.lock()->distance(tgtLoc, srcVec) > ZERO_TOL) {
                        sum += kernel_ptr->evaluate(tgtLoc, srcVec) * srcVals->getScalar(node->coordsContained[i]) *
                            srcWeights->getScalar(node->coordsContained[i]);
                    }
                }
            }
            else {
                for (index_type i = 0; i < node->coordsContained.size(); ++i) {
                    const XyzVector srcVec = _crds.lock()->getVec(node->coordsContained[i]);
                    sum += kernel_ptr->evaluate(tgtLoc, srcVec) * srcVals->getScalar(node->coordsContained[i]) *
                        srcWeights->getScalar(node->coordsContained[i]);
                }
            }
        }
        else {
            for (int i = 0; i < node->kids.size(); ++i) {
                return result + recursiveSum(tgtLoc, srcVals, srcWeights, node->kids[i], tree_depth)
            }
        }
    }
}


}
