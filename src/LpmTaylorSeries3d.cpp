#include "LpmTaylorSeries3d.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>

namespace Lpm {

TaylorCoeffs::TaylorCoeffs(const int maxSeriesOrder) : 
    maxP(maxSeriesOrder) {
    const scalar_type fill_num = std::numeric_limits<scalar_type>::max();
    
    aa.emplace(MultiIndex(0, 0, 0), fill_num);
    
    for (int k1 = 1; k1 <= maxSeriesOrder; ++k1) {
        aa.emplace(MultiIndex(k1, 0, 0), fill_num);
    }
    for (int k2 = 1; k2 <= maxSeriesOrder; ++k2) {
        aa.emplace(MultiIndex(0, k2, 0), fill_num);
    }
        
    for (int k1 = 1; k1 <= maxSeriesOrder; ++k1) {
        for (int k2 = 1; k2 <= maxSeriesOrder - k1; ++k2) {
            aa.emplace(MultiIndex(k1, k2, 0), fill_num);
        }
    }
    
    for (int k3 = 1; k3 <= maxSeriesOrder; ++k3) {
        aa.emplace(MultiIndex(0, 0, k3), fill_num);
        for (int k1 = 1; k1 <= maxSeriesOrder - k3; ++k1) {
            aa.emplace(MultiIndex(k1, 0, k3), fill_num);
        }
        for (int k2 = 1; k2 <= maxSeriesOrder - k3; ++k2) {
            aa.emplace(MultiIndex(0, k2, k3), fill_num);
        }
        
        for (int k1 = 1; k1 <= maxSeriesOrder - k3; ++k1) {
            for (int k2 = 1; k2 <= maxSeriesOrder - k3 - k1; ++k2) {
                aa.emplace(MultiIndex(k1, k2, k3), fill_num);
            }
        }
    }
}

void SphereGreensCoeffs::computeCoeffs(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type param) {
    const scalar_type denom = 1.0 - tgtVec.dotProduct(srcVec);
    aa.at(MultiIndex(0, 0, 0)) = 1.0 / denom;
    
    for (int k1 = 1; k1 <= maxP; ++k1) {
        const scalar_type prev_rec = aa.at(MultiIndex(k1-1, 0, 0));
        const scalar_type newval = tgtVec.x * prev_rec / denom;
        aa.at(MultiIndex(k1, 0, 0)) = newval;
    }
    
    for (int k2 = 1; k2 <= maxP; ++k2) {
        const scalar_type prev_rec = aa.at(MultiIndex(0, k2-1, 0));
        const scalar_type newval = tgtVec.y * prev_rec / denom;
        aa.at(MultiIndex(0, k2, 0)) = newval;
    }
    
    for (int k1 = 1; k1 <= maxP; ++k1) {
        for (int k2 = 1; k2 <= maxP - k1; ++k2) {
            const scalar_type prev_rec = aa.at(MultiIndex(k1, k2-1, 0));
            const scalar_type indexMultiplier = scalar_type(k1 + k2) / k2;
            const scalar_type newval = tgtVec.y * prev_rec * indexMultiplier / denom;
            aa.at(MultiIndex(k1, k2, 0)) = newval;
        }
    }
    
    for (int k3 = 1; k3 <= maxP; ++k3) {
        scalar_type prev_rec = aa.at(MultiIndex(0, 0, k3-1));
        scalar_type newval = tgtVec.z * prev_rec / denom;
        aa.at(MultiIndex(0, 0, k3)) = newval;
        scalar_type indexMultiplier = 1.0;
        
        for (int k1 = 1; k1 <= maxP - k3; ++k1) {
            prev_rec = aa.at(MultiIndex(k1-1, 0, k3));
            indexMultiplier = scalar_type(k1 + k3) / k1;
            newval = tgtVec.x * prev_rec * indexMultiplier / denom;
            aa.at(MultiIndex(k1, 0, k3)) = newval;
        }
        
        for (int k2 = 1; k2 <= maxP - k3; ++k2) {
            prev_rec = aa.at(MultiIndex(0, k2-1, k3));
            indexMultiplier = scalar_type(k2 + k3) / k2;
            newval = tgtVec.y * prev_rec * indexMultiplier / denom;
            aa.at(MultiIndex(0, k2, k3)) = newval;
        }
        
        for (int k1 = 1; k1 <= maxP - k3; ++k1) {
            for (int k2 = 1; k2 <= maxP - k3 - k1; ++k2) {
                prev_rec = aa.at(MultiIndex(k1, k2-1, k3));
                indexMultiplier = scalar_type(k1 + k2 + k3) / k2;
                newval = tgtVec.y * prev_rec * indexMultiplier / denom;
                aa.at(MultiIndex(k1, k2, k3)) = newval;
            }
        }
    }
}
        

std::string TaylorCoeffs::infoString() const {
    std::stringstream ss;
    ss << "TaylorCoeffs info:" << std::endl;
    ss << "\t(multi-index) -- coeff : " << std::endl;
    for (auto& elem : aa) {
        ss << "\t" << elem.first << " -- " << elem.second << std::endl;
    }
    ss << std::endl;
    
    return ss.str();
}



}

