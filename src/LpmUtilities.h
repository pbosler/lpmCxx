#ifndef _LPM_UTILITIES_H
#define _LPM_UTILITIES_H

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include <string>

namespace Lpm {

  scalar_type atan4(const scalar_type y, const scalar_type x);
  
  std::string weight_name(const int ndim);

}
#endif
