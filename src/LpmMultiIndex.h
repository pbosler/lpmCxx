#ifndef _LPM_MULTI_INDEX_H_
#define _LPM_MULTI_INDEX_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include <cmath>
#include <iostream>

namespace Lpm {

inline index_type factorial(const index_type n) {
    return (n == 1 || n == 0) ? 1 : n * factorial(n-1);
}

struct MultiIndex {
    MultiIndex(const unsigned i1 = 0, const unsigned i2 = 0, const unsigned i3 = 0) : k1(i1), k2(i2), k3(i3) {};

    index_type k1;
    index_type k2;
    index_type k3;
    
    inline index_type magnitude() const {return k1 + k2 + k3;}
    
    inline scalar_type vectorPower(const XyzVector& vec) const {
        return std::pow(vec.x, k1) * std::pow(vec.y, k2) * std::pow(vec.z, k3);
    };
    
    inline index_type factorial() const {
        return Lpm::factorial(k1) * Lpm::factorial(k2) * Lpm::factorial(k3);
    }
};

inline bool operator < (const MultiIndex& left, const MultiIndex& right) {
    if (left.k1 != right.k1) 
        return left.k1 < right.k1;
    else {
        if (left.k2 != right.k2)
            return left.k2 < right.k2;
        else
            return left.k3 < right.k3;
    }
}

inline bool operator > (const MultiIndex& left, const MultiIndex& right) {
    return right < left;
}

inline bool operator >= (const MultiIndex& left, const MultiIndex& right) {
    return !(left < right);
}

inline bool operator <= (const MultiIndex& left, const MultiIndex& right) {
    return !(left > right);
}

inline bool operator == (const MultiIndex& left, const MultiIndex& right) {
    return (left.k1 == right.k1 && left.k2 == right.k2) && left.k3 == right.k3;
}

inline bool operator != (const MultiIndex& left, const MultiIndex& right) {
    return !(left == right);
}

inline std::ostream& operator << (std::ostream& os, const MultiIndex& k) {
    os << "(" << k.k1 << ", " << k.k2 << ", " << k.k3 << ")";
    return os;
}



}

#endif
