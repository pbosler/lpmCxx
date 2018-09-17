#include "LpmAosQuadCubicFace.hpp"
#include "LpmGll.hpp"

namespace Lpm {
namespace Aos {

template <int ndim> KidFaceArrays<ndim> QuadCubicFace<ndim>::divide(ParticleSet<ndim>& particles, EdgeSet<ndim>& edges, 
	const index_type myIndex, const index_type faceInsertPoint, 
	const scalar_type radius, const GeometryType geom)  {
	
	KidFaceArrays<ndim> result(12,4,4);
	// parent corners and barycenter
	std::vector<Vec<ndim>> pcorners = this->getCorners(particles);
	std::vector<Vec<ndim>> lcorners = this->getLagCorners(particles);
	Vec<ndim> pctr;
	Vec<ndim> lctr;
	if (geom == PLANAR_GEOMETRY || geom == CARTESIAN_3D_GEOMETRY) {
		pctr = baryCenter(pcorners);
		lctr = baryCenter(lcorners);
	}
	else if (geom == SPHERICAL_SURFACE_GEOMETRY) {
		pctr = sphereBaryCenter(pcorners, radius);
		lctr = sphereBaryCenter(lcorners, radius);
	}
	const index_type ctr_particle = particles.n();
	particles.insert(pctr, lctr, true);
	// connect parent vertices & interiors to child faces
	for (short i=0; i<4; ++i) {
		result.newFaceInteriors[i][i] = this->_interiorInds[i];
		result.newFaceVerts[i][(3*i)%12] = this->_vertInds[(3*i)%12];
		result.newFaceVerts[i][(3*i+1)%12] = this->_vertInds[(3*i+1)%12];
		result.newFaceVerts[i][(3*i+6)%12] = ctr_particle;
		result.newFaceVerts[i][(3*i+11)%12] = this->_vertInds[(3*i+11)%12];
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
			edgeKids[1] = edges.n() + 1;
			edges.divide(parentEdge, particles, radius);
		}
		// connect child edge to new child faces
		if (edges.positiveOrientation(parentEdge, myIndex)) {
			result.newFaceEdges[i][i] = edgeKids[0];
			edges.setLeftFace(edgeKids[0], faceInsertPoint+i);
			result.newFaceVerts[i][(3*i+2)%12] = edges.midpt1(edgeKids[0]);

			result.newFaceEdges[(i+1)%4][i] = edgeKids[1];
			edges.setLeftFace(edgeKids[1], faceInsertPoint+(i+1)%4);
			result.newFaceVerts[i][(3*i+10)%12] = edges.midpt0(edgeKids[1]);
		}
		else {
			result.newFaceEdges[i][i] = edgeKids[1];
			edges.setRightFace(edgeKids[1], faceInsertPoint+i);
			result.newFaceVerts[i][(3*i+2)%12] = edges.midpt0(edgeKids[1]);
			
			result.newFaceEdges[(i+1)%4][i] = edgeKids[0];
			edges.setRightFace(edgeKids[0], faceInsertPoint+(i+1)%4);
			result.newFaceVerts[i][(3*i+10)%12] = edges.midpt1(edgeKids[0]);
		}
		result.newFaceVerts[i][(3*i)%12] = edges.dest(edgeKids[0]);
		result.newFaceVerts[(i+1)%4][(3*i+9)%12] = edges.dest(edgeKids[0]);
	}
	// create new interior edge particles
	std::vector<Vec<ndim>> newpcrds(8);
	std::vector<Vec<ndim>> newlcrds(8);
	const std::vector<index_type> originds = {result.newFaceVerts[0][3], ctr_particle, result.newFaceVerts[1][6], ctr_particle};
	const std::vector<index_type> destinds = {ctr_particle, result.newFaceVerts[3][6], ctr_particle, result.newFaceVerts[3][0]};
	std::vector<Vec<ndim>> origcrds(4);
	std::vector<Vec<ndim>> origlcrds(4);
	std::vector<Vec<ndim>> destcrds(4);
	std::vector<Vec<ndim>> destlcrds(4);
	const index_type particles_insert_point = particles.n();
	for (short i=0; i<4; ++i) {
		origcrds[i] = particles.physCrd(originds[i]);
		origlcrds[i] = particles.lagCrd(originds[i]);
		destcrds[i] = particles.physCrd(destinds[i]);
		destlcrds[i] = particles.lagCrd(destinds[i]);
	}
	if (geom == PLANAR_GEOMETRY || geom == CARTESIAN_3D_GEOMETRY) {
		short j=0;
		for (short i=0; i<4; ++i) {
			newpcrds[j] = pointAlongChord(origcrds[i], destcrds[i], CubicGLL::qp4[1]);
			newlcrds[j++] = pointAlongChord(origlcrds[i], destlcrds[i], CubicGLL::qp4[1]);
			newpcrds[j] = pointAlongChord(origcrds[i], destcrds[i], CubicGLL::qp4[2]);
			newlcrds[j++] = pointAlongChord(origlcrds[i], destlcrds[i], CubicGLL::qp4[2]);
		}
	}
	else if (geom == SPHERICAL_SURFACE_GEOMETRY) {
		short j=0;
		for (short i=0; i<4; ++i) {
			newpcrds[j] = pointAlongCircle(origcrds[i], destcrds[i], CubicGLL::qp4[1], radius);
			newlcrds[j++] = pointAlongCircle(origlcrds[i], destlcrds[i], CubicGLL::qp4[1], radius);
			newpcrds[j] = pointAlongCircle(origcrds[i], destcrds[i], CubicGLL::qp4[2], radius);
			newlcrds[j++] = pointAlongCircle(origlcrds[i], destlcrds[i], CubicGLL::qp4[2], radius);
		}
	}
	for (short i=0; i<8; ++i) {
		particles.insert(newpcrds[i], newlcrds[i], true);
	}
	result.newFaceVerts[0][4] = particles_insert_point;
	result.newFaceVerts[1][11] = particles_insert_point;
	result.newFaceVerts[0][5] = particles_insert_point+1;	
	result.newFaceVerts[1][10] = particles_insert_point+1;
	result.newFaceVerts[3][4] = particles_insert_point+2;
	result.newFaceVerts[2][11] = particles_insert_point+2;
	result.newFaceVerts[3][5] = particles_insert_point+3;
	result.newFaceVerts[2][10] = particles_insert_point+3;
	result.newFaceVerts[1][7] = particles_insert_point+4;
	result.newFaceVerts[2][2] = particles_insert_point+4;
	result.newFaceVerts[1][8] = particles_insert_point+5;
	result.newFaceVerts[2][1] = particles_insert_point+5;
	result.newFaceVerts[0][7] = particles_insert_point+6;
	result.newFaceVerts[3][2] = particles_insert_point+6;
	result.newFaceVerts[0][8] = particles_insert_point+7;
	result.newFaceVerts[3][1] = particles_insert_point+7;
	// debug: make sure all vertices are set
	for (short i=0; i<4; ++i) {
		for (short j=0; j<12; ++j) {
			if (result.newFaceVerts[i][j] < 0) {
				std::ostringstream ss;
				ss << "QuadCubicFace::divide vertex connectivity error, parent " << myIndex
					<< ", child face " << i << ", vertex" << j;
				throw std::runtime_error(ss.str());
			}
		}
	}	
	
	// create new interior edges
	std::vector<std::vector<index_type>> edgeIntInds(4, std::vector<index_type>(2,-1));
	const std::vector<index_type> edgefacelefts = {faceInsertPoint, faceInsertPoint+3, faceInsertPoint+1, faceInsertPoint};
	const std::vector<index_type> edgefacerights = {faceInsertPoint+1, faceInsertPoint+2, faceInsertPoint+2, faceInsertPoint+3};
	edgeIntInds[0] = {particles_insert_point, particles_insert_point+1};
	edgeIntInds[1] = {particles_insert_point+2, particles_insert_point+3};
	edgeIntInds[2] = {particles_insert_point+4, particles_insert_point+5};
	edgeIntInds[3] = {particles_insert_point+6, particles_insert_point+7};
	const index_type edge_insert_point = edges.n();
	for (short i=0; i<4; ++i) {
		edges.insert(originds[i], destinds[i], edgefacelefts[i], edgefacerights[i], edgeIntInds[i]);
	}
	result.newFaceEdges[0][1] = edge_insert_point;
	result.newFaceEdges[1][3] = edge_insert_point;
	result.newFaceEdges[2][3] = edge_insert_point+1;
	result.newFaceEdges[3][1] = edge_insert_point+1;
	result.newFaceEdges[1][2] = edge_insert_point+2;
	result.newFaceEdges[2][0] = edge_insert_point+2;
	result.newFaceEdges[0][2] = edge_insert_point+3;
	result.newFaceEdges[3][0] = edge_insert_point+3;
	// debug: make sure all edges are set
	for (short i=0; i<4; ++i) {
		for (short j=0; j<4; ++j) {
			if (result.newFaceEdges[i][j] < 0) {
				std::ostringstream ss;
				ss << "QuadCubicFace::divide edge connectivity error, parent " << myIndex
					<< ", child face " << i << ", edge " << j;
				throw std::runtime_error(ss.str());
			}
		}
	}
	CubicGLL gll;
	// create new particles for face interiors
	for (short i=0; i<4; ++i) {
		for (short j=0; j<4; ++i) {
			pcorners[j] = particles.physCrd(result.newFaceVerts[i][(3*i)%12]);
			lcorners[j] = particles.lagCrd(result.newFaceVerts[i][(3*i)%12]);
		}
		scalar_type area;
		if (geom == PLANAR_GEOMETRY || geom == CARTESIAN_3D_GEOMETRY) {
			area = polygonArea(pctr, pcorners);
		}
		else if (geom == SPHERICAL_SURFACE_GEOMETRY) {
			area = spherePolygonArea(pctr, pcorners, radius);
		}
		result.kidsFaceArea[i] = area;
		this->_kids[i] = faceInsertPoint+i;
		const std::vector<Vec<ndim>> pintcrds = gll.template quad16interiors<ndim>(pcorners, geom, radius);
		const std::vector<Vec<ndim>> lintcrds = gll.template quad16interiors<ndim>(lcorners, geom, radius);
		particles.move(this->_interiorInds[i], pintcrds[i], lintcrds[i]);
		for (short j=0; j<4; ++j) {
			if (j != i) {
				particles.insert(pintcrds[j], lintcrds[j]);
				result.newFaceInteriors[i][j] = particles.n();
			}
		}
		for (short j=0; j<12; ++j) {
			particles.setWeight(result.newFaceVerts[i][j], area*CubicGLL::quad16edgeqw[j]);
		}
		for (short j=0; j<4; ++j) {
			particles.setWeight(result.newFaceInteriors[i][j], area*CubicGLL::quad16centerqw[j]);
		}
	}
	
	// debug: make sure all interior indices are set
	for (short i=0; i<4; ++i) {
		for (short j=0; j<4; ++j) {
			if (result.newFaceInteriors[i][j] < 0) {
				std::ostringstream ss;
				ss << "QuadCubicFace::divide interior particle connectivity error, parent " << myIndex
					<< ", child " << i << ", index " << j;
				throw std::runtime_error(ss.str());
			}
		}
	}
	this->_area = 0.0;
	return result;
}

template <int ndim> void QuadCubicFace<ndim>::setArea(const scalar_type ar, ParticleSet<ndim>& particles) {
	this->_area = ar;
	for (int j=0; j<12; ++j) {
		particles.setWeight(this->_vertInds[j], ar*CubicGLL::quad16edgeqw[j]);
	}
	for (int j=0; j<4; ++j) {
		particles.setWeight(this->_vertInds[j], ar*CubicGLL::quad16centerqw[j]);
	}
}

template <int ndim> void QuadCubicFace<ndim>::setArea(ParticleSet<ndim>& particles, const GeometryType geom, const scalar_type radius) {
	this->_area = this->computeAreaFromCorners(particles, geom, radius);
	for (int j=0; j<12; ++j) {
		particles.setWeight(this->_vertInds[j], this->_area*CubicGLL::quad16edgeqw[j]);
	}
	for (int j=0; j<4; ++j) {
		particles.setWeight(this->_vertInds[j], this->_area*CubicGLL::quad16centerqw[j]);
	}
}

template <int ndim> std::vector<Vec<ndim>> QuadCubicFace<ndim>::getCorners(const ParticleSet<ndim>& particles) const {
	std::vector<Vec<ndim>> result(4);
	for (short i=0; i<4; ++i) 
		result[i] = particles.physCrd(this->_vertInds[(3*i)%12]);
	return result;
}

template <int ndim> std::vector<Vec<ndim>> QuadCubicFace<ndim>::getLagCorners(const ParticleSet<ndim>& particles) const {
	std::vector<Vec<ndim>> result(4);
	for (short i=0; i<4; ++i) 
		result[i] = particles.lagCrd(this->_vertInds[(3*i)%12]);
	return result;
}



template <int ndim> scalar_type QuadCubicFace<ndim>::computeAreaFromCorners(const ParticleSet<ndim>& particles, const GeometryType geom, 
        	const scalar_type radius) const {
    scalar_type result = 0.0;
    const std::vector<Vec<ndim>> corners = this->getCorners(particles);
    if (geom == PLANAR_GEOMETRY || geom == CARTESIAN_3D_GEOMETRY) {
    	const Vec<ndim> ctr = baryCenter(corners);
    	result = polygonArea(ctr, corners);
    }
    else if (geom == SPHERICAL_SURFACE_GEOMETRY) {
    	const Vec<ndim> ctr = sphereBaryCenter(corners);
    	result = spherePolygonArea(ctr, corners, radius);
    }
    return result;
}

template class QuadCubicFace<2>;
template class QuadCubicFace<3>;
}
}