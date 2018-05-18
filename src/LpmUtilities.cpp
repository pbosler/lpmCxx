#include "LpmUtilities.h"
#include <cmath>

namespace Lpm {

  scalar_type atan4 (const scalar_type y, const scalar_type x) {
    scalar_type result = 0.0;
	if ( x == 0.0 )
	{
		if ( y > 0.0 )
			result = 0.5 * PI;
		else if ( y < 0.0 )
			result = 1.5 * PI;
		else if ( y == 0.0 )
			result = 0.0;
	}
	else if ( y == 0 )
	{
		if ( x > 0.0 )
			result = 0.0;
		else if ( x < 0.0 )
			result = PI;
	}
	else
	{
		scalar_type theta = std::atan2( std::abs(y), std::abs(x) );
		if ( x > 0.0 && y > 0.0 )
			result = theta;
		else if ( x < 0.0 && y > 0.0 )
			result = PI - theta;
		else if ( x < 0.0 && y < 0.0 )
			result = PI + theta;
		else if ( x > 0.0 && y < 0.0 )
			result = 2.0 * PI - theta;
	}
	return result;
  }

}
