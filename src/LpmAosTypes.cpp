#include <iostream>
#include "LpmAosTypes.hpp"
#include "LpmUtilities.h"

#include <cmath>

namespace Lpm {

// scalar_type atan4(const scalar_type y, const scalar_type x) {
//     scalar_type result = 0.0;
// 	if ( x == 0.0 )
// 	{
// 		if ( y > 0.0 )
// 			result = 0.5 * PI;
// 		else if ( y < 0.0 )
// 			result = 1.5 * PI;
// 		else if ( y == 0.0 )
// 			result = 0.0;
// 	}
// 	else if ( y == 0 )
// 	{
// 		if ( x > 0.0 )
// 			result = 0.0;
// 		else if ( x < 0.0 )
// 			result = PI;
// 	}
// 	else
// 	{
// 		scalar_type theta = std::atan2( std::abs(y), std::abs(x) );
// 		if ( x > 0.0 && y > 0.0 )
// 			result = theta;
// 		else if ( x < 0.0 && y > 0.0 )
// 			result = PI - theta;
// 		else if ( x < 0.0 && y < 0.0 )
// 			result = PI + theta;
// 		else if ( x > 0.0 && y < 0.0 )
// 			result = 2.0 * PI - theta;
// 	}
// 	return result;
// }

scalar_type sphereTriArea(const Vec<3>& a, const Vec<3>& b, const Vec<3>& c, const scalar_type radius) {
    const scalar_type s1 = a.sphereDist(b, radius);
    const scalar_type s2 = b.sphereDist(c, radius);
    const scalar_type s3 = c.sphereDist(a, radius);
    const scalar_type halfPerim = 0.5*(s1 + s2 + s3);
    const scalar_type zz = std::tan(0.5*halfPerim) * std::tan(0.5*(halfPerim-s1)) * std::tan(0.5*(halfPerim-s2)) *
        std::tan(0.5*(halfPerim-s3));
    return 4.0 * std::atan(std::sqrt(zz)) * radius*radius;
}

// template <int ndim> std::ostream& operator << (std::ostream& os, const Vec<ndim>& vec) {
//     os << "( "; 
//     for (int i=0; i<ndim; ++i) 
//         os << vec.x[i] << " ";
//     os << " )";
//     return os;
// }

std::ostream& operator << (std::ostream& os, const Vec<1>& vec) {
    os << "(" << vec.x[0] << ")" << std::endl;
}

std::ostream& operator << (std::ostream& os, const Vec<2>& vec) {
    os << "(" << vec.x[0] << ", " << vec.x[1] << ")";
    return os;
}

std::ostream& operator << (std::ostream& os, const Vec<3>& vec) {
    os << "(" << vec.x[0] << ", " << vec.x[1] << ", " << vec.x[2] << ")";
    return os;
}

template class Vec<2>;
template class Vec<3>;

}

