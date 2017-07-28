#include "LpmTaylorSeries3d.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>

namespace Lpm {

TaylorSeries3d::TaylorSeries3d(const int maxSeriesOrder) : 
    maxP(maxSeriesOrder) {
    const scalar_type fill_num = std::numeric_limits<scalar_type>::max();
    
    coeffs.emplace(MultiIndex(0, 0, 0), fill_num);
    moments.emplace(MultiIndex(0,0,0), fill_num);
    
    for (int k1 = 1; k1 <= maxSeriesOrder; ++k1) {
        coeffs.emplace(MultiIndex(k1, 0, 0), fill_num);
        moments.emplace(MultiIndex(k1,0,0), fill_num);
    }
    for (int k2 = 1; k2 <= maxSeriesOrder; ++k2) {
        coeffs.emplace(MultiIndex(0, k2, 0), fill_num);
        moments.emplace(MultiIndex(0,k2,0), fill_num);
    }
        
    for (int k1 = 1; k1 <= maxSeriesOrder; ++k1) {
        for (int k2 = 1; k2 <= maxSeriesOrder - k1; ++k2) {
            coeffs.emplace(MultiIndex(k1, k2, 0), fill_num);
            moments.emplace(MultiIndex(k1,k2,0), fill_num);
        }
    }
    
    for (int k3 = 1; k3 <= maxSeriesOrder; ++k3) {
        coeffs.emplace(MultiIndex(0, 0, k3), fill_num);
        moments.emplace(MultiIndex(0, 0, k3), fill_num);
        for (int k1 = 1; k1 <= maxSeriesOrder - k3; ++k1) {
            coeffs.emplace(MultiIndex(k1, 0, k3), fill_num);
            moments.emplace(MultiIndex(k1,0,k3), fill_num);
        }
        for (int k2 = 1; k2 <= maxSeriesOrder - k3; ++k2) {
            coeffs.emplace(MultiIndex(0, k2, k3), fill_num);
            moments.emplace(MultiIndex(0,k2,k3), fill_num);
        }
        
        for (int k1 = 1; k1 <= maxSeriesOrder - k3; ++k1) {
            for (int k2 = 1; k2 <= maxSeriesOrder - k3 - k1; ++k2) {
                coeffs.emplace(MultiIndex(k1, k2, k3), fill_num);
                moments.emplace(MultiIndex(k1,k2,k3), fill_num);
            }
        }
    }
}

scalar_type TaylorSeries3d::sum() const {
    scalar_type result = 0.0;
    for (auto& elem : coeffs) {
        result += elem.second * moments.at(elem.first);
    }
    return result;
}



void SphereGreensSeries::computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param) {
    const scalar_type denom = 1.0 - tgtVec.dotProduct(srcVec);
    coeffs.at(MultiIndex(0, 0, 0)) = 1.0 / denom;
    
    for (int k1 = 1; k1 <= maxP; ++k1) {
        const scalar_type prev_rec = coeffs.at(MultiIndex(k1-1, 0, 0));
        const scalar_type newval = tgtVec.x * prev_rec / denom;
        coeffs.at(MultiIndex(k1, 0, 0)) = newval;
    }
    
    for (int k2 = 1; k2 <= maxP; ++k2) {
        const scalar_type prev_rec = coeffs.at(MultiIndex(0, k2-1, 0));
        const scalar_type newval = tgtVec.y * prev_rec / denom;
        coeffs.at(MultiIndex(0, k2, 0)) = newval;
    }
    
    for (int k1 = 1; k1 <= maxP; ++k1) {
        for (int k2 = 1; k2 <= maxP - k1; ++k2) {
            const scalar_type prev_rec = coeffs.at(MultiIndex(k1, k2-1, 0));
            const scalar_type indexMultiplier = scalar_type(k1 + k2) / k2;
            const scalar_type newval = tgtVec.y * prev_rec * indexMultiplier / denom;
            coeffs.at(MultiIndex(k1, k2, 0)) = newval;
        }
    }
    
    for (int k3 = 1; k3 <= maxP; ++k3) {
        scalar_type prev_rec = coeffs.at(MultiIndex(0, 0, k3-1));
        scalar_type newval = tgtVec.z * prev_rec / denom;
        coeffs.at(MultiIndex(0, 0, k3)) = newval;
        scalar_type indexMultiplier = 1.0;
        
        for (int k1 = 1; k1 <= maxP - k3; ++k1) {
            prev_rec = coeffs.at(MultiIndex(k1-1, 0, k3));
            indexMultiplier = scalar_type(k1 + k3) / k1;
            newval = tgtVec.x * prev_rec * indexMultiplier / denom;
            coeffs.at(MultiIndex(k1, 0, k3)) = newval;
        }
        
        for (int k2 = 1; k2 <= maxP - k3; ++k2) {
            prev_rec = coeffs.at(MultiIndex(0, k2-1, k3));
            indexMultiplier = scalar_type(k2 + k3) / k2;
            newval = tgtVec.y * prev_rec * indexMultiplier / denom;
            coeffs.at(MultiIndex(0, k2, k3)) = newval;
        }
        
        for (int k1 = 1; k1 <= maxP - k3; ++k1) {
            for (int k2 = 1; k2 <= maxP - k3 - k1; ++k2) {
                prev_rec = coeffs.at(MultiIndex(k1, k2-1, k3));
                indexMultiplier = scalar_type(k1 + k2 + k3) / k2;
                newval = tgtVec.y * prev_rec * indexMultiplier / denom;
                coeffs.at(MultiIndex(k1, k2, k3)) = newval;
            }
        }
    }
}

void TaylorSeries3d::computeMoments(const std::shared_ptr<Coords> crds, const std::vector<index_type>& crdInds, 
            const XyzVector& cntd, const std::shared_ptr<Field> srcVals, const std::shared_ptr<Field> srcWeights) {
    for (auto& elem : moments) {
        const MultiIndex k = elem.first;
        elem.second = 0.0;
        for (index_type i = 0; i < crdInds.size(); ++i) {
            const XyzVector dvec = crds->getVec(crdInds[i]) - cntd;
            const scalar_type strength = srcVals->getScalar(crdInds[i]) * srcWeights->getScalar(crdInds[i]);
            elem.second += k.vectorPower(dvec) * strength;
        }
    }
}

void TaylorSeries3d::computeMoments(const std::shared_ptr<Coords> crds, const std::vector<index_type>& crdInds, 
            const XyzVector& cntd, const std::shared_ptr<Field> srcStrength) {
    for (auto& elem : moments) {
        const MultiIndex k = elem.first;
        elem.second = 0.0;
        for (index_type i = 0; i < crdInds.size(); ++i) {
            const XyzVector dvec = crds->getVec(crdInds[i]) - cntd;
            elem.second += k.vectorPower(dvec) * srcStrength->getScalar(crdInds[i]);
        }
    }
}

void TaylorSeries3d::computeMoments(const std::shared_ptr<Faces> faces, const std::vector<index_type>& crdInds, 
    const XyzVector& cntd, const std::shared_ptr<Field> srcVals) {
    for (auto& elem : moments) {
        const MultiIndex k = elem.first;
        elem.second = 0.0;
        for (index_type i = 0; i < crdInds.size(); ++i) {
            const XyzVector dvec = faces->centroid(crdInds[i]) - cntd;
            elem.second += k.vectorPower(dvec) * srcVals->getScalar(crdInds[i]) * faces->area(crdInds[i]);
        }
    }
}


std::string TaylorSeries3d::infoString() const {
    std::stringstream ss;
    ss << "TaylorSeries info:" << std::endl;
    ss << "\t(multi-index) -- coeff -- moment : " << std::endl;
    for (auto& elem : coeffs) {
        ss << "\tindex " << elem.first << " -- " << elem.second << " -- " << moments.at(elem.first) << std::endl;
    }
    ss << std::endl;
    
    return ss.str();
}

void SecondOrderDelta3dSeries::computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param) {
    const XyzVector dvec = tgtVec - srcVec;
    const scalar_type r = dvec.magnitude();
    const scalar_type expfactor = std::exp(-std::pow(r/param, 3));
    const scalar_type param3 = std::pow(param, 3);
    
    const scalar_type baseval = 3.0 * expfactor / (4.0 * PI * param3);
    
    coeffs.at(MultiIndex(0,0,0)) = baseval;
    coeffs.at(MultiIndex(1,0,0)) = 3.0 * baseval * dvec.x * r / param3;
    coeffs.at(MultiIndex(0,1,0)) = 3.0 * baseval * dvec.y * r / param3;
    coeffs.at(MultiIndex(0,0,1)) = 3.0 * baseval * dvec.z * r / param3;
    coeffs.at(MultiIndex(2,0,0)) = baseval * (-3.0 * dvec.x * dvec.x / (r * param3) - 3.0 * r / param3 + 
        9.0 * dvec.x * dvec.x * r * r / (param3 * param3));
    coeffs.at(MultiIndex(0,2,0)) = baseval * (-3.0 * dvec.y * dvec.y / (r * param3) - 3.0 * r / param3 +
        9.0 * dvec.y * dvec.y * r * r / (param3 * param3));
    coeffs.at(MultiIndex(0,0,2)) = baseval * (-3.0 * dvec.z * dvec.z / (r * param3) - 3.0 * r / param3 +
        9.0 * dvec.z * dvec.z * r * r / (param3 * param3));
    coeffs.at(MultiIndex(1,1,0)) = 9.0 * baseval * dvec.x * dvec.y * r * r / (param3 * param3) - 
        3.0 * baseval * dvec.x * dvec.y / (r * param3);
    coeffs.at(MultiIndex(1,0,1)) = 9.0 * baseval * dvec.x * dvec.z * r * r / (param3 * param3) - 
        3.0 * baseval * dvec.x * dvec.z / (r * param3);
    coeffs.at(MultiIndex(0,1,1)) = 9.0 * baseval * dvec.y * dvec.z * r * r / (param3 * param3) - 
        3.0 * baseval * dvec.y * dvec.z / (r * param3);
}

void SphereDeltaSeries::computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param) {
    const scalar_type bkgrnd = 1.0 / (4.0 * PI);
    const scalar_type dotProd = tgtVec.dotProduct(srcVec);
    const scalar_type distsq = 1.0 - dotProd;
    const scalar_type eps2 = param * param;
    const scalar_type epsdistsq = distsq + eps2;
    const scalar_type baseval = (2.0 * eps2 * dotProd - distsq * distsq ) / (4.0 * PI * epsdistsq * epsdistsq);
    
    coeffs.at(MultiIndex(0,0,0)) = baseval - bkgrnd;
    coeffs.at(MultiIndex(1,0,0)) = 2.0 * tgtVec.x * (eps2 + distsq) / (4.0 * PI * epsdistsq * epsdistsq) +
        tgtVec.x * (-distsq * distsq + 2.0 * eps2 * dotProd) / (2.0 * PI * epsdistsq * epsdistsq * epsdistsq);
    coeffs.at(MultiIndex(0,1,0)) = 2.0 * tgtVec.y * (eps2 + distsq) / (4.0 * PI * epsdistsq * epsdistsq) +
        tgtVec.y * (-distsq * distsq + 2.0 * eps2 * dotProd) / (2.0 * PI * epsdistsq * epsdistsq * epsdistsq);
    coeffs.at(MultiIndex(0,0,1)) = 2.0 * tgtVec.z * (eps2 + distsq) / (4.0 * PI * epsdistsq * epsdistsq) +
        tgtVec.z * (-distsq * distsq + 2.0 * eps2 * dotProd) / (2.0 * PI * epsdistsq * epsdistsq * epsdistsq);
    coeffs.at(MultiIndex(2,0,0)) = (-2.0 * tgtVec.x * tgtVec.x / epsdistsq / epsdistsq + 8.0 * tgtVec.x * tgtVec.x * 
        (eps2 + distsq) / (epsdistsq * epsdistsq * epsdistsq) + 6.0 * tgtVec.x * tgtVec.x * (2.0 * eps2 * dotProd -
        distsq * distsq) / std::pow(epsdistsq, 4)) / (8.0 * PI);
    coeffs.at(MultiIndex(0,2,0)) = (-2.0 * tgtVec.y * tgtVec.y / epsdistsq / epsdistsq + 8.0 * tgtVec.y * tgtVec.y * 
        (eps2 + distsq) / (epsdistsq * epsdistsq * epsdistsq) + 6.0 * tgtVec.y * tgtVec.y * (2.0 * eps2 * dotProd -
        distsq * distsq) / std::pow(epsdistsq, 4)) / (8.0 * PI);
    coeffs.at(MultiIndex(0,0,2)) = (-2.0 * tgtVec.z * tgtVec.z / epsdistsq / epsdistsq + 8.0 * tgtVec.z * tgtVec.z * 
        (eps2 + distsq) / (epsdistsq * epsdistsq * epsdistsq) + 6.0 * tgtVec.z * tgtVec.z * (2.0 * eps2 * dotProd -
        distsq * distsq) / std::pow(epsdistsq, 4)) / (8.0 * PI);
    coeffs.at(MultiIndex(1,1,0)) = - tgtVec.x * tgtVec.y / (2.0 * PI * epsdistsq * epsdistsq) +
        2.0 * tgtVec.y * tgtVec.x * (eps2 + distsq) / (PI * epsdistsq * epsdistsq * epsdistsq) +
        3.0 * tgtVec.x * tgtVec.y * (2.0 * eps2 * dotProd - distsq * distsq) / (2.0 * PI * std::pow(epsdistsq, 4));
    coeffs.at(MultiIndex(1,0,1)) = - tgtVec.x * tgtVec.z / (2.0 * PI * epsdistsq * epsdistsq) +
        2.0 * tgtVec.z * tgtVec.x * (eps2 + distsq) / (PI * epsdistsq * epsdistsq * epsdistsq) +
        3.0 * tgtVec.x * tgtVec.z * (2.0 * eps2 * dotProd - distsq * distsq) / (2.0 * PI * std::pow(epsdistsq, 4));
    coeffs.at(MultiIndex(0,1,1)) = - tgtVec.y * tgtVec.z / (2.0 * PI * epsdistsq * epsdistsq) +
        2.0 * tgtVec.z * tgtVec.y * (eps2 + distsq) / (PI * epsdistsq * epsdistsq * epsdistsq) +
        3.0 * tgtVec.y * tgtVec.z * (2.0 * eps2 * dotProd - distsq * distsq) / (2.0 * PI * std::pow(epsdistsq, 4));
    
}

}

