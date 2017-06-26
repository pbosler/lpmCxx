#ifndef __LPM_EUCLIDEAN_COORDS__
#define __LPM_EUCLIDEAN_COORDS__

#include <iostream>
#include <cmath>
#include <vector>
#include "LpmXyzVector.h"
#include "LpmCoords.h"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"

namespace Lpm {

class EuclideanCoords : public Coords
{
	public :
		EuclideanCoords(const index_type nMax = 0) : Coords(nMax) {};			
		
		scalar_type distance( const index_type indexA, const index_type indexB ) const;
		
		XyzVector midpoint( const index_type indexA, const index_type indexB ) const ;
		
		XyzVector centroid( const std::vector<index_type>& indices ) const;
		
		scalar_type triArea( const index_type indexA, const index_type indexB, const index_type indexC ) const;
};


}
#endif 
