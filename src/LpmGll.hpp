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
//#ifdef USE_APPLE_CLANG
    static constexpr scalar_type sqrt5 = 2.2360679774997896964; // 20-digit precision
//#else
//    static constexpr scalar_type sqrt5 = std::sqrt(5.0);
//#endif
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
    
    /**
        Chooses maximum varying coordinates of a bilinear quadrilateral in R3.  Improves condition number of inverse map. 
    */
    void pickBilinearIJ(int& i, int& j, const std::vector<Aos::Vec<ndim>>& corners) const;
    
    /** Chooses the root in [-1,1], given two roots (r1 and r2) of a quadratic polynomial. */
    scalar_type pickRoot(const scalar_type r1, const scalar_type r2) const;
    
    /** bilinear map of reference element to quadrilateral in R2.
    */
    Aos::Vec<ndim> bilinearMap(const std::vector<Aos::Vec<ndim>>& corners, 
        const scalar_type s1, const scalar_type s2) const;
        
    Aos::Vec<ndim> bilinearMap(const std::vector<Aos::Vec<ndim>>& corners, const Aos::Vec<ndim> refcrds) const {
        return bilinearMap(corners, refcrds[0], refcrds[1]);
    }
    
    /** Jacobian determinant of planar bilinear map */
    scalar_type bilinearMapJacobian(const Aos::Vec<ndim>& refCrds, const std::vector<Aos::Vec<ndim>>& corners) const;
    
    /** Inverts bilinearMap to find reference coordinates (r1,r2) in [-1,1]^2
     Folows C. Hua, 1990, Finite Elem. Anal. Design 7:159--166.
    */
    Aos::Vec<ndim> invertBilinearMap(const std::vector<Aos::Vec<ndim>>& corners, const Aos::Vec<ndim>& queryPt) const;
        
    /** Maps reference element to a bilinear surface in R3, circumscribed by the sphere, then projects that surface onto the sphere.*/
    Aos::Vec<ndim> sphereBilinearMap(const std::vector<Aos::Vec<ndim>>& corners, 
        const scalar_type s1, const scalar_type s2, const scalar_type radius=1.0) const;
        
    inline Aos::Vec<ndim> sphereBilinearMap(const std::vector<Aos::Vec<ndim>>& corners, const Aos::Vec<ndim>& refcrds,
         const scalar_type radius=1.0) const {
         return sphereBilinearMap(corners, refcrds[0], refcrds[1], radius);
    }

    /** Performs cubic interpolation of a scalar using GLL points and basis functions */
    scalar_type quad16interpolation(const Aos::Vec<ndim>& refCrds, const std::vector<scalar_type>& edgeVals, const std::vector<scalar_type>& ctrVals) const;
  };

}
#endif
