#ifndef _LPM_UTILITIES_H
#define _LPM_UTILITIES_H

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include <string>

namespace Lpm {

  /// Inverse tangent with quadrant information, but with output range in [0, 2*pi) instead of (-pi, pi]
  scalar_type atan4(const scalar_type y, const scalar_type x);
  
  /// Returns "length" for 1d, "area" for 2d, "volume" for 3d
  std::string weight_name(const int ndim);
  
  /// Determinant of a 2x2 matrix
  inline scalar_type twoByTwoDeterminant(const scalar_type a, const scalar_type b, const scalar_type c, const scalar_type d) {
    return a*d - b*c;}

   /// Quadratic formula
   void quadraticRoots(scalar_type& r1, scalar_type& r2, const scalar_type a, const scalar_type b, const scalar_type c);
}
#endif
