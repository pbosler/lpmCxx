#include "LpmSakajo2009.h"

namespace Lpm {


SakajoNode::SakajoNode(const Box3d& bbox, SakajoNode* parent, const int max_series_order) : Node(bbox, parent) {
    for (int k1=0; k1<max_series_order; ++k1) {
        for (int k2=0; k2<max_series_order; ++k2) {
            for (int k3=0; k3<max_series_order; ++k3) {
                if (k1 + k2 + k3 <= max_series_order - 1) {
                    const MultiIndex k(k1,k2,k3);
                    coeffs.emplace(k,std::vector<scalar_type>(3, 0.0));
                    moments.emplace(k, std::vector<scalar_type>(3, 0.0));
                }
            }
        }
    }
}


SakajoTree::SakajoTree(const int max_series_order, const int max_tree_depth, const scalar_type sphere_radius, const int prank) : 
Tree(), _maxSeriesOrder(max_series_order), _maxTreeDepth(max_tree_depth), _sphRadius(sphere_radius)
{
    _nnodes = 1;
    log->setProcRank(prank);
    
    const Box3d rbox(sphere_radius);
    std::vector<index_type> rinds;
    
    SakajoNode* rparent = NULL;
    
    _root = std::unique_ptr<SakajoNode>(new SakajoNode(rbox, rparent, max_series_order));
    
    generateTree(dynamic_cast<SakajoNode*>(_root.get()), 0);
}

void SakajoTree::generateTree(SakajoNode* node, const int j) {
    if (j == 3 * _maxSeriesOrder) {
        return;
    }
    else {
        bool split_dim[3];
        for (int i=0; i<3; ++i) {
            split_dim[i] = false;
        }
        split_dim[(j+2)%3] = true;
        
        std::vector<Box3d> kidboxes = node->box.bisectAlongDims(split_dim);
        for (int i=0; i<2; ++i) {
            if (kidboxes[i].intersectsSphere(_sphRadius)) {
                node->kids.push_back(std::unique_ptr<SakajoNode>(new SakajoNode(kidboxes[i], node, _maxSeriesOrder)));
            }
        }
        _nnodes += node->kids.size();
        for (int i=0; i<node->kids.size(); ++i) {
            generateTree(dynamic_cast<SakajoNode*>(node->kids[i].get()), j+1);
        }
    }
}

void SakajoTree::computeCoefficients(const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ) {
    for (index_type i = 0; i < crds->n(); ++i) {
        nodeMoments(dynamic_cast<SakajoNode*>(_root.get()), 0, crds->getVec(i), i, circ->getScalar(i));
    }
}

void SakajoTree::computeCoefficients(const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> vorticity, const std::shared_ptr<Field> area) {
    for (index_type i=0; i<crds->n(); ++i) {
        nodeMoments(dynamic_cast<SakajoNode*>(_root.get()), 0, crds->getVec(i), i, vorticity->getScalar(i) * area->getScalar(i));
    }
}

void SakajoTree::nodeCoeffs(SakajoNode* node, const XyzVector& tgtVec){
    const scalar_type dxy = 1.0 / (_sphRadius*_sphRadius - tgtVec.dotProduct(node->box.centroid()));
    for (int k1=0; k1 < _maxSeriesOrder-1; ++k1) {
        for (int k2=0; k2 < _maxSeriesOrder-1; ++k2) {
            for (int k3=0; k3 < _maxSeriesOrder-1; ++k3) {
                const MultiIndex kkey(k1,k2,k3); 
                const scalar_type num = kkey.magnitude() + 1;
                for (int ii = 0; ii<3; ++ii) {
                    node->coeffs[MultiIndex(k1+1, k2, k3)][ii] = num * dxy * tgtVec.x * node->coeffs[kkey][ii] / (k1 + 1);
                    node->coeffs[MultiIndex(k1, k2+1, k3)][ii] = num * dxy * tgtVec.y * node->coeffs[kkey][ii] / (k2 + 1);
                    node->coeffs[MultiIndex(k1, k2, k3+1)][ii] = num * dxy * tgtVec.z * node->coeffs[kkey][ii] / (k3 + 1);
                }
            }
        }
    }
}


void SakajoTree::nodeMoments(SakajoNode* node, const int k, const XyzVector vecy, const index_type yind, const scalar_type Gamma) {
    if (node->box.containsPoint(vecy)) {
        const XyzVector dvec = vecy - node->box.centroid();
        for (auto& elem : node->moments) {
            const scalar_type vpower = elem.first.vectorPower(dvec);
            elem.second[0] += Gamma * vecy.x * vpower;
            elem.second[1] += Gamma * vecy.y * vpower;
            elem.second[2] += Gamma * vecy.z * vpower;
        }
        if (k == 3 * _maxSeriesOrder) {
            node->coordsContained.push_back(yind);
            return;
        }
        else {
            for (int i=0; i<node->kids.size(); ++i) {
                nodeMoments(dynamic_cast<SakajoNode*>(node->kids[i].get()), k+1, vecy, yind, Gamma);
            }
        }
    }
}

bool SakajoTree::multipoleAcceptance(const SakajoNode* node, const scalar_type meshSize, const scalar_type nuPower, const XyzVector& queryVec) const {
    const scalar_type dist = std::abs(_sphRadius*_sphRadius - queryVec.dotProduct(node->box.centroid()));
    return node->box.radius() <= std::pow(meshSize, nuPower) * dist;
}

XyzVector SakajoTree::biotSavart(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type smooth_param) const {
    XyzVector cp = tgtVec.crossProduct(srcVec);
    const scalar_type denom = (_sphRadius*_sphRadius + smooth_param*smooth_param - tgtVec.dotProduct(srcVec));
    cp.scale(1.0 / denom);
    return cp;
}

void SakajoTree::velocity(XyzVector& vel, SakajoNode* node, const int k, const XyzVector& tgtVec, 
            const scalar_type meshSize, const scalar_type nuPower, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ,
            const scalar_type smooth_param) {
    if (multipoleAcceptance(node, meshSize, nuPower, tgtVec)) {
        nodeCoeffs(node, tgtVec);
        for (auto& elem : node->moments) {
            const std::vector<scalar_type> taylor_coeffs = node->coeffs[elem.first];
            const std::vector<scalar_type> moments = elem.second;
            vel.x += taylor_coeffs[0] * (tgtVec.y * moments[2] - tgtVec.z * moments[1]);
            vel.y += taylor_coeffs[1] * (tgtVec.z * moments[0] - tgtVec.x * moments[2]);
            vel.z += taylor_coeffs[2] * (tgtVec.x * moments[1] - tgtVec.y * moments[0]);
        }
        return;
    }
    else {
        if (k == 3 * _maxSeriesOrder ) {
            for (index_type i=0; i<node->coordsContained.size(); ++i) {
                const XyzVector srcVec = crds->getVec(node->coordsContained[i]);
                const XyzVector kernel = biotSavart(tgtVec, srcVec, smooth_param);
                const scalar_type Gamma = circ->getScalar(node->coordsContained[i]);
                vel.x += Gamma * kernel.x;
                vel.y += Gamma * kernel.y;
                vel.z += Gamma * kernel.z;
            }
            return;
        }
        else {
            for (int i=0; i<node->kids.size(); ++i) {
                velocity(vel, dynamic_cast<SakajoNode*>(node->kids[i].get()), k+1, tgtVec, meshSize, nuPower, crds, circ, smooth_param);
            }
        }
    }
}

}