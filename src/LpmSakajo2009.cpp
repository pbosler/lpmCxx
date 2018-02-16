#include "LpmSakajo2009.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

namespace Lpm {


SakajoNode::SakajoNode(const Box3d& bbox, SakajoNode* parent, const int max_series_order) : Node(bbox, parent) {
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


SakajoTree::SakajoTree(const int max_series_order, const int max_tree_depth, const scalar_type sphere_radius, 
    const scalar_type smooth_param, const int prank) : Tree(), _maxSeriesOrder(max_series_order), 
        _maxTreeDepth(max_tree_depth), _sphRadius(sphere_radius), _smooth(smooth_param)
{
    _nnodes = 1;
    log->setProcRank(prank);
    
    const Box3d rbox(sphere_radius);
    std::vector<index_type> rinds;
    
    SakajoNode* rparent = NULL;
    
    _root = std::unique_ptr<SakajoNode>(new SakajoNode(rbox, rparent, max_series_order));
    
    generateTree(dynamic_cast<SakajoNode*>(_root.get()), 0);
    _depth = computeTreeDepth(_root.get());
}

// Algorithm 1 from Sakajo 2009
void SakajoTree::generateTree(SakajoNode* node, const int j) {
    if (j == 3 * _maxTreeDepth) {
        return;
    }
    else {
        std::vector<bool> split_dim(3, false);
        split_dim[(j+2)%3] = true;
        
        std::vector<Box3d> kidboxes = node->box.bisectAlongDims(split_dim);
        for (int i=0; i<2; ++i) {
            if (kidboxes[i].intersectsSphere()) {
                node->kids.push_back(std::unique_ptr<SakajoNode>(new SakajoNode(kidboxes[i], node, _maxSeriesOrder)));
            }
        }
        _nnodes += node->kids.size();
        for (int i=0; i<node->kids.size(); ++i) {
            generateTree(dynamic_cast<SakajoNode*>(node->kids[i].get()), j+1);
        }
    }
}

void SakajoTree::computeMoments(const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ) {
    for (index_type i = 0; i < crds->n(); ++i) {
        nodeMoments(dynamic_cast<SakajoNode*>(_root.get()), 0, crds->getVec(i), i, circ->getScalar(i));
    }
}

std::vector<MultiIndex> SakajoNode::getKeys() const {
    std::vector<MultiIndex> result;
    for (auto& elem : coeffs) {
        result.push_back(elem.first);
    }
    return result;
}

std::string SakajoNode::coeffString() const {
    std::stringstream ss;
    for (auto& elem : coeffs) {
        ss << elem.first << "  :  " << elem.second << std::endl;
    }
    return ss.str();
}

std::string SakajoNode::infoString() const {
    std::stringstream ss;
    ss << this->Node::infoString();
    ss << "Coeffs" << std::endl;
    ss << coeffString();
    ss << "Moments" << std::endl;
    ss << momentString();
    return ss.str();
}

void SakajoTree::printAll() const {
    printNodeInfo(dynamic_cast<SakajoNode*>(_root.get()));
}

void SakajoTree::printNodeInfo(const SakajoNode* node) const {
    std::cout << node->infoString();
    for (int i=0; i<node->kids.size(); ++i) {
        printNodeInfo(dynamic_cast<SakajoNode*>(node->kids[i].get()));
    }
}

void SakajoTree::writeToVtk(const std::string& filename, const std::string& desc) const {
    std::ofstream fs(filename);
    if (!fs.is_open()) {
        OutputMessage errMsg("cannot open .vtk file", OutputMessage::errorPriority, "Lpm::SakajoTree::writeToVtk");
        log->logMessage(errMsg);
        throw std::ios_base::failure("file write error");  
    }
    fs << "# vtk DataFile Version 2.0" << std::endl;
    fs << desc << std::endl;
    fs << "ASCII" << std::endl;
    fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
    fs << "POINTS " << 8 * _nnodes << " double" << std::endl;
    writeVtkPoints(fs, _root.get());
    fs << "CELLS " << _nnodes << " " << 9 * _nnodes << std::endl;
    index_type vertIndex = 0;
    writeVtkCells(fs, _root.get(), vertIndex);
    fs << "CELL_TYPES " << _nnodes << std::endl;
    writeVtkCellType(fs, _root.get());
    fs << "CELL_DATA " << _nnodes << std::endl;
    fs << "SCALARS tree_level int 1" << std::endl;
    fs << "LOOKUP_TABLE default" << std::endl;
    writeLevelDataToVtk(fs, _root.get());
}

// void SakajoTree::writeMoments(std::ostream& os) const {
//     std::vector<MultiIndex> keys;
//     for (auto& elem : _root->moments) {
//         keys.push_back(elem.first);
//     }    
// }



void SakajoTree::computeMoments(const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> vorticity, const std::shared_ptr<Field> area) {
    for (index_type i=0; i<crds->n(); ++i) {
        nodeMoments(dynamic_cast<SakajoNode*>(_root.get()), 0, crds->getVec(i), i, vorticity->getScalar(i) * area->getScalar(i));
    }
}

// Recursion formulas after equation 11 in Sakajo 2009 
void SakajoNode::calc_coeffs(const XyzVector& tgtVec, const scalar_type sph_radius, const scalar_type delta) {
    const scalar_type denom = square(sph_radius) + square(delta) - tgtVec.dotProduct(box.centroid);
    coeffs.at(MultiIndex(0,0,0)) = 1.0 / denom;
    for (auto& elem : coeffs) {
        const index_type k1 = elem.first.k1;
        const index_type k2 = elem.first.k2;
        const index_type k3 = elem.first.k3;
        auto next_k1 = coeffs.find(MultiIndex(k1+1, k2, k3));
        auto next_k2 = coeffs.find(MultiIndex(k1, k2+1, k3));
        auto next_k3 = coeffs.find(MultiIndex(k1, k2, k3+1));
        
        const scalar_type numer = k1 + k2 + k3 + 1;
      
        if (next_k1 != coeffs.end()) {
            next_k1->second = numer / denom * tgtVec.x * elem.second / (k1 + 1);
            //std::cout << "\tk1: " << elem.second << ", " << MultiIndex(k1+1, k2, k3) << " : " << next_k1->second << std::endl;
        }
        if (next_k2 != coeffs.end()) {
            next_k2->second = numer / denom * tgtVec.y * elem.second / (k2 + 1);
            //std::cout << "\tk2: " << elem.second << ", " << MultiIndex(k1, k2+1, k3) << " : " << next_k1->second << std::endl;
        }
        if (next_k3 != coeffs.end()) {
            next_k3->second = numer / denom * tgtVec.z * elem.second / (k3 + 1);
            //std::cout << "\tk3: " << elem.second << ", " << MultiIndex(k1, k2, k3+1) << " : " << next_k3->second << std::endl;
        }
    }
}


void SakajoNode::calc_moments(const XyzVector& srcY, const scalar_type srcGamma) {
    if (box.containsPoint(srcY)) {
        const XyzVector cntd_to_src = box.centroid - srcY;
        for (auto& elem : moments) {
            const scalar_type vec_power = elem.first.vectorPower(cntd_to_src);
            elem.second[0] += srcGamma * srcY.x * vec_power;
            elem.second[1] += srcGamma * srcY.y * vec_power;
            elem.second[2] += srcGamma * srcY.z * vec_power;
        }
    }
}

std::string SakajoNode::momentString() const {
    std::stringstream ss;
    for (auto& elem : moments) {
        ss << elem.first << " : " << elem.second[0] << "  "  << elem.second[1] << "   " << elem.second[2] << std::endl;
    }
    return ss.str();
}

/* Algorithm 2 from Sakajo 2009, computes equation (13)
    moments[MultiIndex(k)][0] = A_\tau^k
    moments[MultiIndex(k)][1] = B_\tau^k
    moments[MultiIndex(k)][2] = C_\tau^k
*/
void SakajoTree::nodeMoments(SakajoNode* node, const int k, const XyzVector vecy, const index_type yind, const scalar_type Gamma) {
    if (node->box.containsPoint(vecy)) {
        node->calc_moments(vecy, Gamma);
        if (k == 3 * _maxTreeDepth) {
            node->coordsContained.push_back(yind);
        }
        else {
            for (int i=0; i<node->kids.size(); ++i) {
                nodeMoments(dynamic_cast<SakajoNode*>(node->kids[i].get()), k+1, vecy, yind, Gamma);
            }
        }
    }
}

// bool SakajoTree::multipoleAcceptance(const SakajoNode* node, const scalar_type meshSize, const scalar_type nuPower, const XyzVector& queryVec) const {
//     //const scalar_type dist = std::abs(square(_sphRadius) - queryVec.dotProduct(node->box.centroid()));
//     const scalar_type dist = distance(node->box.centroid, queryVec);
//     return node->box.maxRadius <= std::pow(meshSize, nuPower) * dist;
// }

bool SakajoNode::multipoleAcceptance(const XyzVector& queryVec, const scalar_type meshSize, const scalar_type nu) const {
    //return box.maxRadius <= std::pow(meshSize, nu) * distance(box.centroid, queryVec);
    return box.maxRadius / distance(box.centroid, queryVec) < meshSize;
}

XyzVector biotSavart(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type sphRadius, const scalar_type smooth_param) {
    XyzVector cp = tgtVec.crossProduct(srcVec);
    const scalar_type denom = (square(sphRadius) + square(smooth_param) - tgtVec.dotProduct(srcVec)); 
    return cp.scalarMultiply(1.0/denom);
}

scalar_type SakajoTree::greens(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type smooth_param) const {
    return 0.0;
}

/*  Algorithm 3 from Sakajo 2009
*/
void SakajoTree::velocity(XyzVector& vel, SakajoNode* node, const int k, const XyzVector& tgtVec, const index_type tgtInd,
            const scalar_type meshSize, const scalar_type nuPower, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ) {
    if (node->multipoleAcceptance(tgtVec, meshSize, nuPower)) {
        node->calc_coeffs(tgtVec, _sphRadius, _smooth);
        for (auto& elem : node->coeffs) {
            const std::vector<scalar_type> mm = node->moments.at(elem.first);
            vel.x += elem.second * (tgtVec.y * mm[2] - tgtVec.z * mm[1]);
            vel.y += elem.second * (tgtVec.z * mm[0] - tgtVec.x * mm[2]);
            vel.z += elem.second * (tgtVec.x * mm[1] - tgtVec.y * mm[0]);
        }
    }
    else {
        if (k == 3 * _maxTreeDepth ) {
            for (index_type i=0; i<node->coordsContained.size(); ++i) {
                const index_type srcInd = node->coordsContained[i];
                if (tgtInd != srcInd) {
                    const XyzVector srcVec = crds->getVec(srcInd);
                    const XyzVector kernel = biotSavart(tgtVec, srcVec, 1.0, _smooth);
                    const scalar_type Gamma = circ->getScalar(srcInd);
                    vel += kernel.scalarMultiply(Gamma);
                }
            }
            return;
        }
        else {
            for (int i=0; i<node->kids.size(); ++i) {
                velocity(vel, dynamic_cast<SakajoNode*>(node->kids[i].get()), k+1, tgtVec, tgtInd, meshSize, nuPower, crds, circ);
            }
        }
    }
}

void SakajoTree::computeVelocity(std::shared_ptr<Field> outputVelocity, const std::shared_ptr<SphericalCoords> crds, 
    const std::shared_ptr<Field> circ, const scalar_type meshSize, const scalar_type nuPower) {
    for (index_type i=0; i<crds->n(); ++i) {
        XyzVector vel(0.0, 0.0, 0.0);
        const XyzVector tgtVec = crds->getVec(i);
        velocity(vel, dynamic_cast<SakajoNode*>(_root.get()), 0, tgtVec, i, meshSize, nuPower, crds, circ); 
        outputVelocity->replace(i, vel.scalarMultiply(-1.0/(4.0 * PI * _sphRadius)));
    }
}

// void SakajoTree::streamFn(scalar_type& psi, SakajoNode* node, const int k, const XyzVector& tgtVec, 
//             const scalar_type meshSize, const scalar_type nuPower, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ) {
//     if (multipoleAcceptance(node, meshSize, nuPower, tgtVec)) {
//         nodeCoeffs(node, tgtVec);
//         for (auto& elem : node->moments) {
//             const scalar_type taylor_coeff = node->coeffs[elem.first];
//             for (int ii=0; ii<3; ++ii) {
//                 psi += taylor_coeff * elem.second[ii];
//             }
//             psi /= (- 4.0 * PI * _sphRadius);
//             return;
//         }
//     }
//     else {
//         if (k==3 * _maxSeriesOrder) {
//             for (index_type i = 0; i < node->coordsContained.size(); ++i) {
//                 const XyzVector srcVec = crds->getVec(node->coordsContained[i]);
//                 const scalar_type kernel = greens(tgtVec, srcVec, _smooth);
//                 const scalar_type Gamma = circ->getScalar(node->coordsContained[i]);
//                 psi += Gamma * kernel;
//             }
//             return;
//         }
//         else {
//             for (int i=0; i<node->kids.size(); ++i) {
//                 streamFn(psi, dynamic_cast<SakajoNode*>(node->kids[i].get()), k+1, tgtVec, meshSize, nuPower, crds, circ);
//             }
//         }
//     }
// }

}