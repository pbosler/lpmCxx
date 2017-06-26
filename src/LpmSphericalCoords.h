#ifndef __LPM_SphereCoords__
#define __LPM_SphereCoords__

#include <iostream>
#include <cmath>
#include <vector>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmCoords.h"

namespace Lpm {

class SphericalCoords : public Coords {
	public :
		SphericalCoords( const index_type nMax = 0,	const scalar_type radius = 1.0) : 
			Coords(nMax), _radius(radius) {};
		
		scalar_type distance( const index_type indexA, const index_type indexB ) const;
		
		XyzVector midpoint( const index_type indexA, const index_type indexB ) const; 
		
		XyzVector centroid( const std::vector<index_type>& indices ) const;
		
		scalar_type triArea( const index_type indexA, const index_type indexB, const index_type indexC ) const;
		
		scalar_type Latitude( const index_type index ) const;
		
		scalar_type Longitude( const index_type index ) const;
		
		inline scalar_type radius() const {return _radius;}
		
	protected : 
		scalar_type _radius;
};


}
#endif

