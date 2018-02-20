#include "LpmTreeSum.h"
#include <cmath>

namespace Lpm { 

SumNode::SumNode(const Box3d& bbox, SumNode* parent, const std::vector<index_type>& crd_indices, 
    const int max_series_order) : Node(bbox, parent, crd_indices) {
    for (int k1=0; k1<max_series_order; ++k1) {
        for (int k2=0; k2<max_series_order; ++k2) {
            for (int k3=0; k3<max_series_order; ++k3) {
                if (k1 + k2 + k3 < max_series_order) {
                    const MultiIndex k(k1,k2,k3);
                    coeffs.emplace(k, 0.0);
                    moments.emplace(k, std::vector<scalar_type>(3, 0.0));
                }
            }
        }
    }
}

bool SumNode::multipoleAcceptance(const XyzVector& tgtVec, const scalar_type tol) const {
    return box.maxRadius / distance(box.centroid, tgtVec) < tol;
}

void TreeSum::vecSum(XyzVector& vel, const XyzVector& tgtVec, SumNode* node, const Coords* crds, 
    const scalar_type mp_tol, const std::shared_ptr<Field> strength) {
    if (node->multipoleAcceptance(tgtVec, mp_tol)) {
        node->calc_coeffs(tgtVec, _smooth);
        for (auto& elem : node->coeffs) {
            const std::vector<scalar_type> mm = node->moments.at(elem.first);
            vel.x += elem.second * (tgtVec.y * mm[2] - tgtVec.z * mm[1]);
            vel.y += elem.second * (tgtVec.z * mm[0] - tgtVec.x * mm[2]);
            vel.z += elem.second * (tgtVec.x * mm[1] - tgtVec.y * mm[0]);
        }
    }
    else {
        if (node->isLeaf()) {
            for (index_type j=0; j<node->coordsContained.size(); ++j) {
                const XyzVector srcVec = crds->getVec(node->coordsContained[j]);
                const scalar_type str = strength->getScalar(node->coordsContained[j]);
                XyzVector kernel = vKern->evaluate(tgtVec, srcVec);
                kernel.scale(str);
                vel += kernel;
            }
        }
        else {
            for (int i=0; i<node->kids.size(); ++i) {
                vecSum(vel, tgtVec, dynamic_cast<SumNode*>(node->kids[i].get()), crds, mp_tol, strength);
            }
        }
    }
}

void SumNode::calc_moments(const XyzVector& srcVec, const scalar_type srcStrength) {
    if (box.containsPoint(srcVec)) {
        const XyzVector cntd_to_src = box.centroid - srcVec;
        for (auto & elem : moments ) {
            const scalar_type vecPower = elem.first.vectorPower(cntd_to_src);
            elem.second[0] += srcStrength * srcVec.x * vecPower;
            elem.second[1] += srcStrength * srcVec.y * vecPower;
            elem.second[2] += srcStrength * srcVec.z * vecPower;
        }
    }
}

void SumNode::calc_coeffs(const XyzVector& tgtVec, const scalar_type smoother) {
    const scalar_type denom = 1.0 + square(smoother) - tgtVec.dotProduct(box.centroid);
    coeffs.at(MultiIndex(0,0,0)) = 1.0 / denom;
    for (auto& elem : coeffs) {
        const index_type k1 = elem.first.k1;
        const index_type k2 = elem.first.k2;
        const index_type k3 = elem.first.k3;
        auto nextK1 = coeffs.find(MultiIndex(k1+1, k2, k3));
        auto nextK2 = coeffs.find(MultiIndex(k1, k2+1, k3));
        auto nextK3 = coeffs.find(MultiIndex(k1, k2, k3+1));
        
        const scalar_type numer = k1 + k2 + k3 + 1;
        if (nextK1 != coeffs.end()) {
            nextK1->second = numer / denom * tgtVec.x * elem.second / (k1+1);
        }
        if (nextK2 != coeffs.end()) {
            nextK2->second = numer / denom * tgtVec.y * elem.second / (k2 + 1);
        }
        if (nextK3 != coeffs.end()) {
            nextK3->second = numer / denom * tgtVec.z * elem.second / (k3 + 1);
        }
    }
}

void TreeSum::computeMoments(const std::shared_ptr<Field> strength) {
    std::shared_ptr<Coords> crds = _crds.lock();
    for (index_type i=0; i<crds->n(); ++i) {
        const XyzVector& srcVec = crds->getVec(i);
        const scalar_type str = strength->getScalar(i);
        nodeMoments(dynamic_cast<SumNode*>(_root.get()), srcVec, str);
    }
}

void TreeSum::nodeMoments(SumNode* node, const XyzVector& srcVec, const scalar_type srcStrength) {
    if (node->box.containsPoint(srcVec)) {
        node->calc_moments(srcVec, srcStrength);
        for (int i=0; i<node->kids.size(); ++i) {
            nodeMoments(dynamic_cast<SumNode*>(node->kids[i].get()), srcVec, srcStrength);
        }
    }
}



TreeSum::TreeSum(const std::shared_ptr<Coords> crds, const int max_series_order, const int max_particles_per_box, 
    const scalar_type smoother, const bool do_shrink) :
    Tree(crds), _maxSeriesOrder(max_series_order), _maxParticlesPerNode(max_particles_per_box), _smooth(smoother), 
    _shrink(do_shrink), sKern(NULL), vKern(NULL)
{
    _nnodes = 1;
    const Box3d rbox = _root->box;
    const std::vector<index_type>* inds = &_root->coordsContained;
    _root = std::unique_ptr<SumNode>(new SumNode(rbox, NULL, *inds, _maxSeriesOrder));
    
    generateTree(dynamic_cast<SumNode*>(_root.get()), _crds.lock().get());
}

void TreeSum::generateTree(SumNode* node, const Coords* crds) {
    if (node->coordsContained.size() > _maxParticlesPerNode) {
        std::vector<Box3d> kidboxes = node->box.bisectAll();
        for (int i=0; i<kidboxes.size(); ++i) {
            std::vector<index_type> kidcoords;
            kidcoords.reserve(node->coordsContained.size());  
            for (index_type j=0; j<node->coordsContained.size(); ++j) {
                if (kidboxes[i].containsPoint(crds->getVec(node->coordsContained[j]))) {
                    kidcoords.push_back(node->coordsContained[j]);
                }
            }
            if (!kidcoords.empty()) {
                kidcoords.shrink_to_fit();
                node->kids.push_back(std::unique_ptr<SumNode>(new SumNode(kidboxes[i], node, kidcoords, _maxSeriesOrder)));
            } 
        }
        _nnodes += node->kids.size();
        for (int i=0; i<node->kids.size(); ++i) {
            if (_shrink) {
                shrinkBox(node->kids[i].get());
            }
            generateTree(dynamic_cast<SumNode*>(node->kids[i].get()), crds);
        }
    }
}

}

