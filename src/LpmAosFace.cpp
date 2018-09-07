#include <iostream>

#include <sstream>
#include "LpmAosFace.hpp"
#include "LpmAosEdge.hpp"

namespace Lpm {
namespace Aos {

template <int ndim> std::string KidFaceArrays<ndim>::infoString() const {
	std::ostringstream ss;
	ss << "new face data:"<< std::endl;
	for (short i=0; i<4; ++i) {
		ss << "Face " << i << ":" << std::endl;
		ss << "vertices/edge particles: ";
		for (short j=0; j<newFaceVerts.size(); ++j)
			ss << newFaceVerts[i][j] << " ";
		ss << std::endl << "edges: ";
		for (short j=0; j<newFaceEdges.size(); ++j) 
			ss << newFaceEdges[i][j] << " ";
		ss << std::endl << "interior particles: ";
		for (short j=0; j<newFaceInteriors.size(); ++j)
			ss << newFaceInteriors[i][j] << " ";
		ss << std::endl << "vertex weights: ";
		for (short j=0; j<newVertWeights.size(); ++j) 
			ss << newVertWeights[i][j] << " ";
		ss << std::endl << "interior weights: ";
		for (short j=0; j<newInteriorWeights.size(); ++j) 
			ss << newInteriorWeights[i][j] << " ";
		ss << std::endl;
	}
	return ss.str();
}

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
		result.newInteriorWeights[i][0] = ar;
	}
	this->_area = 0.0;
	particles.setWeight(this->_interiorInds[0], 0.0);
	return result;
}

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
			edgeKids[1] = edges.n() + 1;
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
		result.newFaceVerts[i][(i+1)%3] = edges.dest(edgeKids[1]);
		result.newFaceVerts[(i+1)%3][i] = edges.dest(edgeKids[1]);
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
	edges.insert(result.newFaceVerts[4][0], result.newFaceVerts[4][1], faceInsertPoint+3, faceInsertPoint+2);
	result.newFaceEdges[2][0] = edge_insert_point;
	result.newFaceEdges[3][0] = edge_insert_point;
	
	edges.insert(result.newFaceVerts[4][1], result.newFaceVerts[4][2], faceInsertPoint+3, faceInsertPoint);
	result.newFaceEdges[0][1] = edge_insert_point+1;
	result.newFaceEdges[3][1] = edge_insert_point+1;
	
	edges.insert(result.newFaceVerts[4][2], result.newFaceVerts[4][0], faceInsertPoint+3, faceInsertPoint+1);
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
		result.newInteriorWeights[i][0] = ar;
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
	particles.move(this->_interiorInds[0], pcenter, lcenter);
	particles.setWeight(this->_interiorInds[0], ar);
	this->_area = 0.0;
	return result;
}

template <int ndim> KidFaceArrays<ndim> QuadCubicFace<ndim>::divide(ParticleSet<ndim>& particles, EdgeSet<ndim>& edges, 
	const index_type myIndex, const index_type faceInsertPoint, 
	const scalar_type radius, const GeometryType geom)  {
	KidFaceArrays<ndim> result;
	return result;
}

template <int ndim> Vec<ndim> Face<ndim>::physBarycenter(const ParticleSet<ndim>& particles) const {
    Vec<ndim> result;
    for (int i=0; i<_vertInds.size(); ++i) {
        result += particles.physCrd(_vertInds[i]);
    }
    result.scaleInPlace(1.0/_vertInds.size());
    return result;
}

template <int ndim> Vec<ndim> Face<ndim>::lagBarycenter(const ParticleSet<ndim>& particles) const {
    Vec<ndim> result;
    for (int i=0; i<_vertInds.size(); ++i) {
        result += particles.lagCrd(_vertInds[i]);
    }
    result.scaleInPlace(1.0/_vertInds.size());
    return result;
}

template <int ndim> Vec<3> Face<ndim>::physSphBarycenter(const ParticleSet<ndim>& particles, 
    const scalar_type radius) const {
    Vec<3> result;
    for (int i=0; i<_vertInds.size(); ++i) {
        result += particles.physCrd(_vertInds[i]);
    }
    result.scaleInPlace(1.0/_vertInds.size());
    result.scaleInPlace(radius);
    return result;
}

template <int ndim> Vec<3> Face<ndim>::lagSphBarycenter(const ParticleSet<ndim>& particles, 
    const scalar_type radius) const {
    Vec<3> result;
    for (int i=0; i<_vertInds.size(); ++i) {    
        result += particles.lagCrd(_vertInds[i]);
    }
    result.scaleInPlace(1.0/_vertInds.size());
    result.scaleInPlace(radius);
    return result;
}

template <int ndim> scalar_type Face<ndim>::scalarIntegral(const std::string field_name, const ParticleSet<ndim>& particles) const {
	scalar_type result = 0.0;
	for (int i=0; i<_vertInds.size(); ++i) {
		result += particles.scalarVal(_vertInds[i], field_name)*particles.weight(_vertInds[i]);
	}
	for (int i=0; i<_interiorInds.size(); ++i) {
		result += particles.scalarVal(_interiorInds[i], field_name)*particles.weight(_interiorInds[i]);
	}
	return result;
}

template <int ndim> std::string Face<ndim>::infoString() const {
    std::ostringstream ss;
    ss << "Face info:" << std::endl;
    ss << "\tvertInds = ";
    for (int i=0; i<_vertInds.size(); ++i) 
        ss << _vertInds[i] << " ";
    ss << std::endl;
    ss << "\tedgeInds = ";
    for (int i=0; i<_edgeInds.size(); ++i)
        ss << _edgeInds[i] << " ";
    ss << std::endl;
    ss << "\tinterior indices = " ;
    for (int i=0; i<_interiorInds.size(); ++i) 
        ss << _interiorInds[i] << " ";
    ss << std::endl;
    ss << "\tparent face = " << _parent << std::endl;
    ss << "\tkids = ";
    for (int i=0; i<4; ++i)
        ss << _kids[i] << " ";
    ss << std::endl;
    ss << "\tarea = " << _area << std::endl;
    return ss.str();
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

template <int ndim> std::vector<Vec<ndim>> TriFace<ndim>::getCorners(const ParticleSet<ndim>& particles) const {
	std::vector<Vec<ndim>> result(3);
	for (short i=0; i<3; ++i) 
		result[i] = particles.physCrd(this->_vertInds[i]);
	return result;
}

template <int ndim> std::vector<Vec<ndim>> QuadCubicFace<ndim>::getCorners(const ParticleSet<ndim>& particles) const {
	std::vector<Vec<ndim>> result(4);
	for (short i=0; i<4; ++i) 
		result[i] = particles.physCrd(this->_vertInds[(3*i)%12]);
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

template class Face<2>;
template class Face<3>;
template class QuadFace<2>;
template class QuadFace<3>;
template class TriFace<2>;
template class TriFace<3>;
template class QuadCubicFace<2>;
template class QuadCubicFace<3>;
template struct KidFaceArrays<2>;
template struct KidFaceArrays<3>;
}
}
