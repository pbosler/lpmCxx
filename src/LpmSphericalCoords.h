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
			Coords(nMax, true), _radius(radius) {_geometry = SPHERICAL_SURFACE_GEOMETRY;}
		
		scalar_type distance( const index_type indexA, const index_type indexB ) const;
		scalar_type distance(const XyzVector& v0, const index_type ind) const;
		scalar_type distance(const XyzVector& v0, const XyzVector& v1) const;
		
		XyzVector midpoint( const index_type indexA, const index_type indexB ) const; 
		
		XyzVector centroid( const std::vector<index_type>& indices ) const;
		
		scalar_type triArea( const index_type indexA, const index_type indexB, const index_type indexC ) const;
		
		scalar_type triArea( const XyzVector& v0, const index_type indexA, const index_type indexB ) const;
		
		scalar_type latitude( const index_type ind ) const;
		
		scalar_type longitude( const index_type ind ) const;
		
		inline scalar_type radius() const {return _radius;}
		
		void initRandom(const bool useTimeSeed = false, const scalar_type domainRadius = 1.0);
		
	protected : 
		scalar_type _radius;
};


}
#endif

