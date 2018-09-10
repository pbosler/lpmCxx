#include "LpmAosTriFace.hpp"

namespace Lpm {
namespace Aos {

template <int ndim> KidFaceArrays<ndim> TriFace<ndim>::divide(ParticleSet<ndim>& particles, EdgeSet<ndim>& edges, 
	const index_type myIndex, const index_type faceInsertPoint, 
	const scalar_type radius, const GeometryType geom)  {
	KidFaceArrays<ndim> result(3,3,1);
	// parent vertices
	for (short i=0; i<3; ++i) {
		result.newFaceVerts[i][i] = this->_vertInds[i];
	}
	// loop over parent edges
	for (short i=0; i<3; ++i) {
		const index_type parentEdge = this->_edgeInds[i];
		std::array<index_type, 2> edgeKids;
		// check if edge has already been divided (by an adjacent panel)
		if (edges.isDivided(parentEdge)) {
			edgeKids = edges.kids(parentEdge);
		}
		else {
			edgeKids[0] = edges.n();
			edgeKids[1] = edges.n()+1;
			edges.divide(parentEdge, particles, radius);
		}
		// connect child edges to new faces
		if (edges.positiveOrientation(parentEdge, myIndex)) {
			result.newFaceEdges[i][i] = edgeKids[0];
			edges.setLeftFace(edgeKids[0], faceInsertPoint+i);
			
			result.newFaceEdges[(i+1)%3][i] = edgeKids[1];
			edges.setLeftFace(edgeKids[1], faceInsertPoint+i+1);
		}
		else {
			result.newFaceEdges[i][i] = edgeKids[1];
			edges.setRightFace(edgeKids[1], faceInsertPoint+i);
			
			result.newFaceEdges[(i+1)%3][i] = edgeKids[0];
			edges.setRightFace(edgeKids[0], faceInsertPoint+i+1);
		}
		result.newFaceVerts[i][(i+1)%3] = edges.dest(edgeKids[0]);
		result.newFaceVerts[(i+1)%3][i] = edges.dest(edgeKids[0]);
	}
	// TriFace only: construct child 3
	result.newFaceVerts[3][0] = result.newFaceVerts[2][1];
	result.newFaceVerts[3][1] = result.newFaceVerts[2][0];
	result.newFaceVerts[3][2] = result.newFaceVerts[1][0];
	// debug: make sure all vertices are set
	for (short i=0; i<4; ++i) {
		for (short j=0; j<3; ++j) {
			if (result.newFaceVerts[i][j] < 0) {
				std::ostringstream ss;
				ss << "TriFace::divide vertex connectivity error, parent " << myIndex
					<< ", child face " << i << ", vertex" << j;
				throw std::runtime_error(ss.str());
			}
		}
	}
	
	// create new interior edges
	const index_type edge_insert_point = edges.n();
	edges.insert(result.newFaceVerts[3][0], result.newFaceVerts[3][1], faceInsertPoint+3, faceInsertPoint+2);
	result.newFaceEdges[2][0] = edge_insert_point;
	result.newFaceEdges[3][0] = edge_insert_point;
	
	edges.insert(result.newFaceVerts[3][1], result.newFaceVerts[3][2], faceInsertPoint+3, faceInsertPoint);
	result.newFaceEdges[0][1] = edge_insert_point+1;
	result.newFaceEdges[3][1] = edge_insert_point+1;
	
	edges.insert(result.newFaceVerts[3][2], result.newFaceVerts[3][0], faceInsertPoint+3, faceInsertPoint+1);
	result.newFaceEdges[1][2] = edge_insert_point+2;
	result.newFaceEdges[3][2] = edge_insert_point+2;
	// debug: make sure all edges are set
	for (short i=0; i<4; ++i) {
		for (short j=0; j<3; ++j) {
			if (result.newFaceEdges[i][j] < 0) {
				std::ostringstream ss;
				ss << "TriFace::divide edge connectivity error, parent " << myIndex
					<< ", child face " << i << ", edge " << j;
				throw std::runtime_error(ss.str());
			}
		}
	}
	
	// create new particles for face centers
	const index_type particle_insert_point = particles.n();
	std::vector<Vec<ndim>> pcorners(3);
	std::vector<Vec<ndim>> lcorners(3);
	Vec<ndim> pcenter;
	Vec<ndim> lcenter;
	scalar_type ar;
	for (short i=0; i<3; ++i) {
		for (short j=0; j<3; ++j) {
			pcorners[j] = particles.physCrd(result.newFaceVerts[i][j]);
			lcorners[j] = particles.lagCrd(result.newFaceVerts[i][j]);
		}
		if (geom == PLANAR_GEOMETRY || geom == CARTESIAN_3D_GEOMETRY) {
			pcenter = baryCenter(pcorners);
			lcenter = baryCenter(lcorners);
			ar = polygonArea(pcenter, pcorners);
		}
		else if (geom == SPHERICAL_SURFACE_GEOMETRY) {
			pcenter = sphereBaryCenter(pcorners, radius);
			lcenter = sphereBaryCenter(lcorners, radius);
			ar = spherePolygonArea(pcenter, pcorners, radius);
		}
		particles.insert(pcenter, lcenter, ar);
		result.newFaceInteriors[i][0] = particle_insert_point+i;
		result.kidsFaceArea[i] = ar;
		this->_kids[i] = faceInsertPoint+i;
	}
	// TriFace only: move parent center particle to child 3 center particle
	for (short j=0; j<3; ++j) {
		pcorners[j] = particles.physCrd(result.newFaceVerts[3][j]);
		lcorners[j] = particles.lagCrd(result.newFaceVerts[3][j]);
	}
	if (geom == PLANAR_GEOMETRY || geom == CARTESIAN_3D_GEOMETRY) {
		pcenter = baryCenter(pcorners);
		lcenter = baryCenter(lcorners);
		ar = polygonArea(pcenter, pcorners);
	}
	else if (geom == SPHERICAL_SURFACE_GEOMETRY) {
		pcenter = sphereBaryCenter(pcorners, radius);
		lcenter = sphereBaryCenter(lcorners, radius);
		ar = spherePolygonArea(pcenter, pcorners, radius);
	}
	result.newFaceInteriors[3][0] = this->_interiorInds[0];
	particles.move(this->_interiorInds[0], pcenter, lcenter);
	particles.setWeight(this->_interiorInds[0], ar);
	this->_kids[3] = faceInsertPoint+3;
	result.kidsFaceArea[3] = ar;
	this->_area = 0.0;
	return result;
}

template <int ndim> scalar_type TriFace<ndim>::computeAreaFromCorners(const ParticleSet<ndim>& particles, const GeometryType geom, 
        	const scalar_type radius) const {
    scalar_type result = 0.0;
    const Vec<ndim> ctr = particles.physCrd(this->_interiorInds[0]);
    const std::vector<Vec<ndim>> corners = this->getCorners(particles);
    if (geom == PLANAR_GEOMETRY || geom == CARTESIAN_3D_GEOMETRY) {
    	result = polygonArea(ctr, corners);
    }
    else if (geom == SPHERICAL_SURFACE_GEOMETRY) {
    	result = spherePolygonArea(ctr, corners, radius);
    }
    return result;
}

template <int ndim> std::vector<Vec<ndim>> TriFace<ndim>::getCorners(const ParticleSet<ndim>& particles) const {
	std::vector<Vec<ndim>> result(3);
	for (short i=0; i<3; ++i) 
		result[i] = particles.physCrd(this->_vertInds[i]);
	return result;
}

template <int ndim> std::vector<Vec<ndim>> TriFace<ndim>::getLagCorners(const ParticleSet<ndim>& particles) const {
	std::vector<Vec<ndim>> result(3);
	for (short i=0; i<3; ++i) 
		result[i] = particles.lagCrd(this->_vertInds[i]);
	return result;
}


template class TriFace<2>;
template class TriFace<3>;

}
}