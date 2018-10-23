#ifndef LPM_GLL_HPP
#define LPM_GLL_HPP

#include "LpmTypeDefs.h"
#include "LpmConfig.h"
#include "LpmAosTypes.hpp"
#include "LpmUtilities.h"
#include <cmath>
#include <array>

namespace Lpm {

template <int ndim=2>  struct CubicGLL {
    static constexpr scalar_type sqrt5 = std::sqrt(5.0);
    static constexpr scalar_type oosqrt5 = 1.0/sqrt5;
    static constexpr scalar_type oo6 = 1.0/6.0;
    static constexpr scalar_type qp0 = -1.0;
    static constexpr scalar_type qp1 = -oosqrt5;
    static constexpr scalar_type qp2 = oosqrt5;
    static constexpr scalar_type qp3 = 1.0;
    static constexpr scalar_type qw0 = oo6;
    static constexpr scalar_type qw1 = 1.0-oo6;
    static constexpr scalar_type qw2 = 1.0-oo6;
    static constexpr scalar_type qw3 = oo6;
        
    inline scalar_type basis0(const scalar_type s) const {return 0.125*(-1.0 + s + 5.0*s*s - 5.0*s*s*s);}
    inline scalar_type basis1(const scalar_type s) const {return 0.625*(s-1.0)*(s+1.0)*(sqrt5*s - 1.0);}
    inline scalar_type basis2(const scalar_type s) const {return -0.125*sqrt5*(s-1.0)*(s+1.0)*(sqrt5 + 5.0*s);}
    inline scalar_type basis3(const scalar_type s) const {return 0.125*(-1 - s + 5*s*s +5*s*s*s);}
    
    scalar_type qp4(const int id) const;
    scalar_type qw4(const int id) const;
    
    Aos::Vec<ndim> quad16edgeqp(const int id) const;
    Aos::Vec<ndim> quad16centerqp(const int id) const;
    scalar_type quad16edgeqw(const int id) const;
    scalar_type quad16centerqw(const int id) const;
    
    void pickBilinearIJ(int& i, int& j, const std::vector<Aos::Vec<ndim>>& corners) const;
    scalar_type pickRoot(const scalar_type r1, const scalar_type r2) const;
    
    Aos::Vec<ndim> invertBilinearMap(const std::vector<Aos::Vec<ndim>>& corners, const Aos::Vec<ndim>& queryPt) const;
    
    Aos::Vec<ndim> bilinearMap(const std::vector<Aos::Vec<ndim>>& corners, 
        const scalar_type s1, const scalar_type s2) const;
        
    scalar_type bilinearMapJacobian(const Aos::Vec<ndim>& refCrds, const std::vector<Aos::Vec<ndim>>& corners) const;

    Aos::Vec<ndim> sphereBilinearMap(const std::vector<Aos::Vec<ndim>>& corners, 
        const scalar_type s1, const scalar_type s2, const scalar_type radius=1.0) const;
    
    std::vector<Aos::Vec<ndim>> quad16interiors(const std::vector<Aos::Vec<ndim>>& corners, 
        const GeometryType geom, const scalar_type radius = 1.0) const;
        
    scalar_type quad16interpolation(const Aos::Vec<ndim>& refCrds, const std::vector<scalar_type>& edgeVals, const std::vector<scalar_type>& ctrVals) const;
  };

}
#endif
