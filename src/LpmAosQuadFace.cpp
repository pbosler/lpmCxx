#include "LpmAosQuadFace.hpp"

namespace Lpm {
namespace Aos {

template <int ndim> KidFaceArrays<ndim> QuadFace<ndim>::divide(ParticleSet<ndim>& particles, EdgeSet<ndim>& edges, 
	const index_type myIndex, const index_type faceInsertPoint, 
	const scalar_type radius, const GeometryType geom)  {
	KidFaceArrays<ndim> result(4,4,1);
	// parent vertices and center particle
	// QuadFace only: change parent center particle to a vertex of each child panel
	for (short i=0; i<4; ++i) {
		result.newFaceVerts[i][i] = this->_vertInds[i];
		result.newFaceVerts[i][(i+2)%4] = this->_interiorInds[0];
	}
	// loop over parent edges
	for (short i=0; i<4; ++i) {
		const index_type parentEdge = this->_edgeInds[i];
		std::array<index_type,2> edgeKids;
		// Check if edge has already been divided (by an adjacent face)
		if (edges.isDivided(parentEdge)) {
			edgeKids = edges.kids(parentEdge);
		}
		else {
			edgeKids[0] = edges.n();
			edgeKids[1] = edges.n()+1;
			edges.divide(parentEdge, particles, radius);
		}
		// connect child edges to new child faces
		if (edges.positiveOrientation(parentEdge, myIndex)) {
			result.newFaceEdges[i][i] = edgeKids[0];
			edges.setLeftFace(edgeKids[0], faceInsertPoint+i);
			
			result.newFaceEdges[(i+1)%4][i] = edgeKids[1];
			edges.setLeftFace(edgeKids[1], faceInsertPoint+i+1); 
		}
		else{
			result.newFaceEdges[i][i] = edgeKids[1];
			edges.setRightFace(edgeKids[1], faceInsertPoint+i);
			
			result.newFaceEdges[(i+1)%4][i] = edgeKids[0];
			edges.setRightFace(edgeKids[0], faceInsertPoint+i+1);
		}
		result.newFaceVerts[i][(i+1)%4] = edges.dest(edgeKids[1]);
		result.newFaceVerts[(i+1)%4][i] = edges.dest(edgeKids[1]);
	}
	// debug: make sure all vertices are set
	for (short i=0; i<4; ++i) {
		for (short j=0; j<4; ++j) {
			if (result.newFaceVerts[i][j] < 0) {
				std::ostringstream ss;
				ss << "QuadFace::divide vertex connectivity error, parent " << myIndex
					<< ", child face " << i << ", vertex" << j;
				throw std::runtime_error(ss.str());
			}
		}
	}	
	// create new interior edges
	const index_type edge_insert_point = edges.n();
	edges.insert(result.newFaceVerts[0][1], result.newFaceVerts[0][2], faceInsertPoint, faceInsertPoint+1);
	result.newFaceEdges[0][1] = edge_insert_point;
	result.newFaceEdges[1][3] = edge_insert_point;
	
	edges.insert(result.newFaceVerts[3][1], result.newFaceVerts[3][2], faceInsertPoint+3, faceInsertPoint+2);
	result.newFaceEdges[2][3] = edge_insert_point+1;
	result.newFaceEdges[3][1] = edge_insert_point+1;
	
	edges.insert(result.newFaceVerts[2][1], result.newFaceVerts[2][0], faceInsertPoint+1, faceInsertPoint+2);
	result.newFaceEdges[1][2] = edge_insert_point+2;
	result.newFaceEdges[2][0] = edge_insert_point+2;
	
	edges.insert(result.newFaceVerts[0][2], result.newFaceVerts[0][3], faceInsertPoint, faceInsertPoint+3);
	result.newFaceEdges[0][2] = edge_insert_point+3;
	result.newFaceEdges[3][0] = edge_insert_point+3;
	
	// debug: make sure all edges are set
	for (short i=0; i<4; ++i) {
		for (short j=0; j<4; ++j) {
			if (result.newFaceEdges[i][j] < 0) {
				std::ostringstream ss;
				ss << "QuadFace::divide edge connectivity error, parent " << myIndex
					<< ", child face " << i << ", edge " << j;
				throw std::runtime_error(ss.str());
			}
		}
	}
	
	// create new particles for face centers
	const index_type particle_insert_point = particles.n();
	std::vector<Vec<ndim>> pcorners(4);
	std::vector<Vec<ndim>> lcorners(4);
	for (short i=0; i<4; ++i) {
		for (short j=0; j<4; ++j) {
			pcorners[j] = particles.physCrd(result.newFaceVerts[i][j]);
			lcorners[j] = particles.lagCrd(result.newFaceVerts[i][j]);
		}
		Vec<ndim> pcenter;
		Vec<ndim> lcenter;
		scalar_type ar;
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
	this->_area = 0.0;
	particles.setWeight(this->_interiorInds[0], 0.0);
	return result;
}

template <int ndim> scalar_type QuadFace<ndim>::computeAreaFromCorners(const ParticleSet<ndim>& particles, const GeometryType geom, 
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

template <int ndim> std::vector<Vec<ndim>> QuadFace<ndim>::getCorners(const ParticleSet<ndim>& particles) const {
	std::vector<Vec<ndim>> result(4);
	for (short i=0; i<4; ++i) 
		result[i] = particles.physCrd(this->_vertInds[i]);
	return result;
}

template <int ndim> std::vector<Vec<ndim>> QuadFace<ndim>::getLagCorners(const ParticleSet<ndim>& particles) const {
	std::vector<Vec<ndim>> result(4);
	for (short i=0; i<4; ++i) 
		result[i] = particles.lagCrd(this->_vertInds[i]);
	return result;
}

template class QuadFace<2>;
template class QuadFace<3>;

}
}