#ifndef LPM_GLL_HPP
#define LPM_GLL_HPP

#include "LpmTypeDefs.h"
#include "LpmConfig.h"
#include "LpmAosTypes.hpp"
#include "LpmUtilities.h"
#include <cmath>

namespace Lpm {

  struct CubicGLL {
    static constexpr scalar_type sqrt5 = std::sqrt(5.0);
    static constexpr scalar_type oosqrt5 = 1.0/sqrt5;
    static constexpr scalar_type oo6 = 1.0/6.0;
    static constexpr scalar_type qp4[4] = {-1.0, -oosqrt5, oosqrt5, 1.0};
    static constexpr scalar_type qw4[4] = {oo6, 1.0-oo6, 1.0-oo6, oo6};
    static constexpr scalar_type quad16edgeqp[24] = 
        {qp4[0], qp4[3], qp4[0], qp4[2], qp4[0], qp4[1], qp4[0], qp4[0], // edge0, (x,y) coord pairs
         qp4[1], qp4[0], qp4[2], qp4[0], qp4[3], qp4[0], // edge1
         qp4[3], qp4[1], qp4[3], qp4[2], qp4[3], qp4[3], // edge2
         qp4[2], qp4[3], qp4[1], qp4[3]}; // edge3
    static constexpr scalar_type quad16edgeqw[12] = {
         qw4[0]*qw4[3], qw4[0]*qw4[2], qw4[0]*qw4[1], qw4[0]*qw4[0], // edge0
         qw4[1]*qw4[0], qw4[2]*qw4[0], qw4[3]*qw4[0], // edge1
         qw4[3]*qw4[1], qw4[3]*qw4[2], qw4[3]*qw4[3], // edge2
         qw4[2]*qw4[3], qw4[1]*qw4[3]};
    static constexpr scalar_type quad16centerqp[8] = {
         qp4[1], qp4[2], qp4[1], qp4[1], qp4[2], qp4[1], qp4[2], qp4[2]};
    static constexpr scalar_type quad16centerqw[4] = {
        qw4[1]*qw4[2], qw4[1]*qw4[1], qw4[2]*qw4[1], qw4[2]*qw4[2]};
    
    inline scalar_type basis0(const scalar_type s) const {return 0.125*(-1.0 + s + 5.0*s*s - 5.0*s*s*s);}
    inline scalar_type basis1(const scalar_type s) const {return 0.625*(s-1.0)*(s+1.0)*(sqrt5*s - 1.0);}
    inline scalar_type basis2(const scalar_type s) const {return -0.125*sqrt5*(s-1.0)*(s+1.0)*(sqrt5 + 5.0*s);}
    inline scalar_type basis3(const scalar_type s) const {return 0.125*(-1 - s + 5*s*s +5*s*s*s);}
  };

}
#endif
