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
		EuclideanCoords(const index_type nMax = 0, const GeometryType gkind = PLANAR_GEOMETRY) : 
		    Coords(nMax, (gkind == PLANAR_GEOMETRY ? false : true)) {_geometry = gkind;}
		    
		void initRandom(const bool useTimeSeed = false, const scalar_type domainRadius = 1.0);
		
		scalar_type distance( const index_type indexA, const index_type indexB ) const;
		scalar_type distance(const XyzVector& vec, const index_type ind) const;
		scalar_type distance(const XyzVector& v0, const XyzVector& v1) const;
		
		XyzVector midpoint( const index_type indexA, const index_type indexB ) const ;
		
		XyzVector centroid( const std::vector<index_type>& indices ) const;
		
		scalar_type triArea( const index_type indexA, const index_type indexB, const index_type indexC ) const;
		
		scalar_type triArea( const XyzVector& v0, const index_type indexA, const index_type indexB ) const;
};


}
#endif 
