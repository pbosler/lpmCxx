#include "LpmSakajo2009.h"

namespace Lpm {


SakajoNode::SakajoNode(const Box3d& bbox, SakajoNode* parent, const int max_series_order) : Node(bbox, parent) {
    for (int k1=0; k1<max_series_order; ++k1) {
        for (int k2=0; k2<max_series_order; ++k2) {
            for (int k3=0; k3<max_series_order; ++k3) {
                if (k1 + k2 + k3 <= max_series_order - 1) {
                    const MultiIndex k(k1,k2,k3);
                    coeffs.emplace(k,std::vector<scalar_type>(3, 0.0));
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

void SakajoTree::computeCoefficients(const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> field) {
    for (index_type i = 0; i < crds->n(); ++i) {
        nodeCoeffs(dynamic_cast<SakajoNode*>(_root.get()), 0, crds->getVec(i), i, field->getScalar(i));
    }
}


void SakajoTree::nodeCoeffs(SakajoNode* node, const int k, const XyzVector vecy, const index_type yind, const scalar_type Gamma) {
    if (node->box.containsPoint(vecy)) {
        const XyzVector dvec = vecy - node->box.centroid();
        for (auto& elem : node->coeffs) {
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
                nodeCoeffs(dynamic_cast<SakajoNode*>(node->kids[i].get()), k+1, vecy, yind, Gamma);
            }
        }
    }
}

}