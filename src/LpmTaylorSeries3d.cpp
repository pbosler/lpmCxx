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


}

