#include "LpmTreeSum.h"
#include "LpmMultiIndex.h"
#include "LpmTaylorSeries3d.h"
#include "LpmOutputMessage.h"
#include "LpmBox3d.h"
#include <cmath>
#include <sstream>
#include <limits>
#include <exception>
#include <mpi.h>

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

void SumNode::computeMoments(const std::shared_ptr<Faces> faces, const std::shared_ptr<Field> srcVals) {
    series->computeMoments(faces, coordsContained, box.centroid(), srcVals);
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

TreeSum::TreeSum(std::shared_ptr<MeshedParticles> pmesh, const std::shared_ptr<ScalarKernel> kernel, 
    const std::string tgtFieldName, const std::string srcFieldName, const scalar_type maxAspectRatio, 
    const int maxSeriesOrder, const scalar_type seriesParam, const int prank, const scalar_type farFieldParam) :
    Tree(), _maxP(maxSeriesOrder), _kernel(kernel), _seriesParam(seriesParam), _nuParam(farFieldParam), 
    _mesh(pmesh->meshPtr()), _vertTgt(pmesh->getVertexFieldPtr(tgtFieldName)), 
    _faceTgt(pmesh->getFaceFieldPtr(tgtFieldName)), _faceSrc(pmesh->getFaceFieldPtr(srcFieldName)) 
{
    log->setProcRank(prank);
    
    std::shared_ptr<Coords> crds = _mesh.lock()->getPhysCoords();
    Box3d rbox(crds->minX() - BOX_PADDING_FACTOR * crds->minX(), 
      crds->maxX() + BOX_PADDING_FACTOR * crds->maxX(),
      crds->minY() - BOX_PADDING_FACTOR * crds->minY(), 
      crds->maxY() + BOX_PADDING_FACTOR * crds->maxY(), 
      crds->minZ() - BOX_PADDING_FACTOR * crds->minZ(), 
      crds->maxZ() + BOX_PADDING_FACTOR * crds->maxZ());
      
    std::vector<index_type> face_inds;
    std::shared_ptr<Faces> faces = _mesh.lock()->getFaces();
    face_inds.reserve(faces->nLeaves());
    for (index_type i = 0; i < faces->n(); ++i) {
        if (faces->isLeaf(i)) {
            face_inds.push_back(i);
        }
    }
    
    Node* rparent = NULL;
    
    _root = std::unique_ptr<Node>(new SumNode(rbox, rparent, face_inds, _maxP, kernel.get()));
    
    std::stringstream ss;
    ss << "Tree constructed for MeshedParticles" << std::endl;
    ss << infoString();
    ss << "\tmax source value = " << _faceSrc.lock()->maxScalarVal() << std::endl;
    ss << "\tmin source value = " << _faceSrc.lock()->minScalarVal() << std::endl;
    ss << "\tmesh surface area = " << _mesh.lock()->surfaceArea() << std::endl;
    OutputMessage traceMsg(ss.str(), OutputMessage::tracePriority, "TreeSum::TreeSum");
    log->logMessage(traceMsg);
}


void TreeSum::buildTree(const index_type maxCoordsPerNode) {
    generateTree(_root.get(), maxCoordsPerNode);
    _depth = computeTreeDepth(_root.get());
}

void TreeSum::buildTreeFromMesh(const index_type maxCoordsPerNode) {
    generateTreeFromMesh(_root.get(), maxCoordsPerNode);
    _depth = computeTreeDepth(_root.get());
}

void TreeSum::generateTreeFromMesh(Node* node, const index_type maxCoordsPerNode) {
    if (node->coordsContained.size() <= maxCoordsPerNode) {
        return;
    }
    else {
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
        
        std::shared_ptr<Faces> faces = _mesh.lock()->getFaces();
        std::vector<Box3d> kidboxes = node->box.bisectAlongDims(splitDims);
        for (int i = 0; i < kidboxes.size(); ++i) {
            std::vector<index_type> kidcoords;
            kidcoords.reserve(node->coordsContained.size());
            for (index_type j = 0; j < node->coordsContained.size(); ++j) {
                if (kidboxes[i].containsPoint(faces->centroid(node->coordsContained[j]))) {
                    kidcoords.push_back(node->coordsContained[j]);
                }
            }
            if (!kidcoords.empty()) {
                kidcoords.shrink_to_fit();
                node->kids.push_back(std::unique_ptr<Node>(new SumNode(kidboxes[i], node, kidcoords, _maxP, _kernel.lock().get())));
            }
        }
        if (node->kids.empty()) {
            OutputMessage errMsg("All kids are empty, this shouldn't happen.", OutputMessage::errorPriority, "Treee::generateTreeFromMesh");
            log->logMessage(errMsg);
            throw std::runtime_error("Tree error encountered");
        }
        else {
            _nnodes += node->kids.size();
            for (int i = 0; i < node->kids.size(); ++i) {
                shrinkBox(node->kids[i].get(), faces);  
                generateTreeFromMesh(node->kids[i].get(), maxCoordsPerNode);
            }
        }
    }
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
               recursiveSum(sum, tgtLoc, srcVals, srcWeights, dynamic_cast<SumNode*>(node->kids[i].get()), tree_depth);
            }
        }
    }
}

scalar_type TreeSum::recursiveSumMeshVertices(const index_type& tgtInd, std::shared_ptr<Field> srcVals, 
    SumNode* node, const index_type& tree_depth) 
{
    scalar_type result = 0.0;
    const XyzVector tgtLoc = _mesh.lock()->getPhysCoords()->getVec(tgtInd);
    if (node->isFar(tgtLoc, tree_depth, _nuParam)) {
        if (!node->momentsReady) {
            node->computeMoments(_mesh.lock()->getFaces(), srcVals);
        }
        node->computeCoeffs(tgtLoc, _seriesParam);
        result += node->seriesSum();
        //std::cout << "series sum = " << node->seriesSum() << std::endl;
    }
    else {
        if (node->isLeaf()) {
            std::shared_ptr<ScalarKernel> kernel = _kernel.lock();
            std::shared_ptr<Faces> faces = _mesh.lock()->getFaces();
            scalar_type dval = 0.0;
            std::stringstream ss;
            //ss << "(kernel, srcval, area) = [";
            for (index_type i = 0; i < node->coordsContained.size(); ++i) {
                const XyzVector srcLoc = faces->centroid(node->coordsContained[i]);
                scalar_type kval = kernel->evaluate(tgtLoc, srcLoc);
                scalar_type sval = srcVals->getScalar(node->coordsContained[i]);
                scalar_type aval = faces->area(node->coordsContained[i]);
                dval += kval * sval * aval;
                //ss << "(" << kval << ", " << sval << ", " << aval << ")  ";
            }
            //ss << "]" << std::endl;;
            //ss << ", local sum = " << dval << std::endl;
            //std::cout << ss.str();
            
            result += dval;
        }
        else {
            for (int i = 0; i < node->kids.size(); ++i) {
                result += recursiveSumMeshVertices(tgtInd, srcVals, dynamic_cast<SumNode*>(node->kids[i].get()), tree_depth);
            }
        }
    }
    return result;
}

scalar_type TreeSum::recursiveSumMeshFaces(const index_type& tgtInd, std::shared_ptr<Faces> faces, std::shared_ptr<Field> srcVals,
    SumNode* node, const index_type& tree_depth)
{
    scalar_type result = 0.0;
    const XyzVector tgtLoc = faces->centroid(tgtInd);
    if (node->isFar(tgtLoc, tree_depth, _nuParam)) {
        if (!node->momentsReady) {
            node->computeMoments(faces, srcVals);
        }
        node->computeCoeffs(tgtLoc, _seriesParam);
        result += node->seriesSum();
    }
    else {
        if (node->isLeaf()) {
            std::shared_ptr<ScalarKernel> kernel = _kernel.lock();
            if (kernel->isSingular()) {
                for (index_type i = 0; i < node->coordsContained.size(); ++i) {
                    if ( tgtInd != node->coordsContained[i]) {
                        const XyzVector srcLoc = faces->centroid(node->coordsContained[i]);
                        result += kernel->evaluate(tgtLoc, srcLoc) * srcVals->getScalar(node->coordsContained[i]) * 
                            faces->area(node->coordsContained[i]);
                    }
                }
            }
            else {
                for (index_type i = 0; i < node->coordsContained.size(); ++i) {
                    const XyzVector srcLoc = faces->centroid(node->coordsContained[i]);
                    result += kernel->evaluate(tgtLoc, srcLoc) * srcVals->getScalar(node->coordsContained[i]) * 
                        faces->area(node->coordsContained[i]);
                }
            }
        }
        else {
            for (int i = 0; i < node->kids.size(); ++i) {
                result += recursiveSumMeshFaces(tgtInd, faces, srcVals, dynamic_cast<SumNode*>(node->kids[i].get()), tree_depth);
            }
        }
    }
    return result;
}


void TreeSum::meshSolve(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) {
    const int myRank = mpiVerts.getRank();
    const int commSize = mpiVerts.getSize();
    
    std::shared_ptr<Faces> face_ptr = _mesh.lock()->getFaces();
    std::shared_ptr<Field> src_ptr = _faceSrc.lock();
    std::shared_ptr<Field> vert_tgt_ptr = _vertTgt.lock();
    std::shared_ptr<Field> face_tgt_ptr = _faceTgt.lock();

    const int tree_depth = depth();
    
    for (index_type i = mpiVerts.startIndex(myRank); i <= mpiVerts.endIndex(myRank); ++i) {
        vert_tgt_ptr->replace(i, recursiveSumMeshVertices(i, src_ptr, dynamic_cast<SumNode*>(_root.get()), tree_depth));
    }
    
    for (index_type i = mpiFaces.startIndex(myRank); i <= mpiFaces.endIndex(myRank); ++i) {
        if (face_ptr->isLeaf(i)) {
            face_tgt_ptr->replace(i, recursiveSumMeshFaces(i, face_ptr, src_ptr, dynamic_cast<SumNode*>(_root.get()), tree_depth));
        }
    }
}

void TreeSum::meshBroadcast(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const {
    for (int i = 0; i < mpiVerts.getSize(); ++i) {
        scalar_type* bufStart0 = _vertTgt.lock()->getPtrToData();
        scalar_type* bufStart1 = _faceTgt.lock()->getPtrToData();
        MPI_Bcast(bufStart0 + mpiVerts.startIndex(i), mpiVerts.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
        MPI_Bcast(bufStart1 + mpiFaces.startIndex(i), mpiFaces.msgSize(i), MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
}

}
