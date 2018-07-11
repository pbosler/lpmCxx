#include "LpmAosFaceSet.hpp"
#include <exception>
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace Lpm {

template <int ndim> void FaceSet<ndim>::insert(const ind_vec& intrs, const ind_vec& verts, const ind_vec& edges, 
            const index_type pt, const scalar_type ar) {
    if (_faces.size() + 1 >= _nMax) {
        throw std::out_of_range("FaceSet::insert _nMax exceeded.");
    }
    _faces.push_back(_factory->createFace(intrs, verts, edges, pt, ar));
    _nActive += 1;
}

template <int ndim> scalar_type FaceSet<ndim>::minArea() const {
    scalar_type result = std::numeric_limits<scalar_type>::max();
    for (index_type i=0; i<_faces.size(); ++i) {
        if (_faces[i]->area() < result)
            result = _faces[i]->area();
    }
    return result;
}

template <int ndim> scalar_type FaceSet<ndim>::maxArea() const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_faces.size(); ++i) {
        if (_faces[i]->area() > result)
            result = _faces[i]->area();
    }
    return result;
}

template <int ndim> scalar_type FaceSet<ndim>::maxLeafArea() const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_faces.size(); ++i) {
        if (_faces[i]->isLeaf()) {
            if (_faces[i]->area() > result) 
                result = _faces[i]->area();
        }
    }
    return result;
}

template <int ndim> std::string FaceSet<ndim>::infoString(const bool printAll) const {
    std::ostringstream ss;
    ss << "FaceSet info:" << std::endl;
    ss << "\tgeom = " << geometryString(_geom) << std::endl;
    ss << "\tnMax = " << _nMax << std::endl;
    ss << "\tsize = " << _faces.size() << std::endl;
    ss << "\tnActive = " << _nActive << std::endl;
    if (printAll) {
        for (index_type i=0; i<_faces.size(); ++i) {
            ss << _faces[i]->infoString();
        }
    }
    return ss.str();
}


template <int ndim> void FaceSet<ndim>::divide(const index_type ind, ParticleSet<ndim>& particles, EdgeSet<ndim>& edges) {
    if (_faces.size() + 4 >= _nMax) {
        throw std::out_of_range("FaceSet::divide nMax exceeded.");
    }
    
    const std::vector<index_type> pfaceVerts = _faces[ind]->vertices();
    const std::vector<index_type> pfaceEdges = _faces[ind]->edges();
    const std::vector<index_type> pInteriors = _faces[ind]->interiors();
    if (dynamic_cast<TriFace<ndim>*>(_faces[0].get())) {
        std::vector<std::vector<index_type>> newFaceVerts(4, std::vector<index_type>(3));
        std::vector<std::vector<index_type>> newFaceEdges(4, std::vector<index_type>(3));
        std::array<index_type, 4> newFaceInds;
        std::array<index_type, 3> newEdgeInds;
        for (int i=0; i<3; ++i) {
            newEdgeInds[i] = edges.n() + i;
        }        
        for (int i=0; i<4; ++i) {
            newFaceInds[i] = _faces.size()+i;
        }
        
        Vec<ndim> bctr;
        Vec<ndim> lagbctr;
        if (_geom == SPHERICAL_SURFACE_GEOMETRY) {
            bctr = _faces[ind]->physSphBarycenter(particles, _radius);
            lagbctr = _faces[ind]->lagSphBarycenter(particles, _radius);
        }
        else {
            bctr = _faces[ind]->physBarycenter(particles);
            lagbctr = _faces[ind]->lagBarycenter(particles);
        }
        
        for (int i=0; i<3; ++i) {
            newFaceVerts[i][i] = pfaceVerts[i];
        }
        
        // loop over parent edges
        for (int i=0; i<3; ++i) {
            const index_type parentEdge = pfaceEdges[i];
            std::array<index_type,2> childedges;
            if (edges.isDivided(ind)) {
                childedges = edges.kids(parentEdge);
            }
            else {
                childedges[0] = edges.n();
                childedges[1] = edges.n() + 1;
                edges.divide(parentEdge, particles);
            }
            
            newFaceVerts[i][(i+1)%3] = edges.dest(childedges[0]);
            newFaceVerts[(i+1)%3][i] = edges.dest(childedges[0]);
            newFaceVerts[3][(i+2)%3] = edges.dest(childedges[0]);
            
            if (edges.positiveOrientation(parentEdge, ind)) {
                newFaceEdges[i][i] = childedges[0];
                edges.setLeftFace(childedges[0], newFaceInds[i]);
                
                newFaceEdges[(i+1)%3][i] = childedges[1];
                edges.setLeftFace(childedges[1], newFaceInds[(i+1)%3]);
            }
            else {
                newFaceEdges[i][i] = childedges[1];
                edges.setRightFace(childedges[1], newFaceInds[i]);
                
                newFaceEdges[(i+1)%3][i] = childedges[0];
                edges.setRightFace(childedges[0], newFaceInds[(i+1)%3]);
            }
        }
        
        
        edges.insert(newFaceVerts[2][1], newFaceVerts[2][0], newFaceInds[3], newFaceInds[2]);
        newFaceEdges[3][0] = newEdgeInds[0];
        newFaceEdges[2][0] = newEdgeInds[0];
        
        edges.insert(newFaceVerts[0][2], newFaceVerts[0][1], newFaceInds[3], newFaceInds[0]);
        newFaceEdges[3][1] = newEdgeInds[1];
        newFaceEdges[0][1] = newEdgeInds[1];
        
        edges.insert(newFaceVerts[1][2], newFaceVerts[1][0], newFaceInds[3], newFaceInds[1]);
        newFaceEdges[3][2] = newEdgeInds[2];
        newFaceEdges[1][2] = newEdgeInds[2];
        
        // create new center particles
        std::vector<Vec<ndim>> triCorners(3);
        std::vector<Vec<ndim>> lagtriCorners(3);
        std::vector<Vec<ndim>> newFaceCenters(4);
        std::vector<Vec<ndim>> newFaceLagCenters(3);
        for (int i=0; i<4; ++i) {
            for (int j=0; j<3; ++j) {
                triCorners[j] = particles.physCrd(newFaceVerts[i][j]);
                lagtriCorners[j] = particles.lagCrd(newFaceVerts[i][j]);
            }   
            if (_geom == SPHERICAL_SURFACE_GEOMETRY) {
                newFaceCenters[i] = sphereBaryCenter(triCorners, _radius);
                newFaceLagCenters[i] = sphereBaryCenter(lagtriCorners, _radius);
            }
            else {
                newFaceCenters[i] = baryCenter(triCorners);
                newFaceLagCenters[i] = baryCenter(lagtriCorners);
            }    
        }
        // create child faces
        const index_type faceInsert = _faces.size();
        for (int i=0; i<3; ++i) {
            particles.insert(newFaceCenters[i], newFaceLagCenters[i]);
            this->insert(std::vector<index_type>(particles.n()-1), newFaceVerts[i], newFaceEdges[i], ind);
            this->_faces[ind]->_kids[i] = faceInsert  + i;
        }
        // reuse parent face's center particle for child face 3
        this->insert(pInteriors, newFaceVerts[3], newFaceEdges[3], ind);
        particles.move(pInteriors[0], newFaceCenters[3], newFaceLagCenters[3]);
        
        for (index_type i=faceInsert; _faces.size(); ++i) {
            this->_faces[i]->setArea(_geom, particles, _radius);
        }
    }
    else if (dynamic_cast<QuadFace<ndim>*>(_faces[0].get())) {
        std::vector<std::vector<index_type>> newFaceVerts(4, std::vector<index_type>(4));
        std::vector<std::vector<index_type>> newFaceEdges(4, std::vector<index_type>(4));
        
        // connect parent vertices to child faces
        // convert center particle to vertex particle
        for (int i=0; i<4; ++i) {
            newFaceVerts[i][i] = pfaceVerts[i];
            newFaceVerts[i][(i+2)%4] = pInteriors[0];
        }
        
        // loop over parent edges, divide edges if required
        for (int i=0; i<4; ++i) {
            const index_type parentEdge = pfaceEdges[i];
            std::array<index_type,2> childedges;
            if (edges.isDivided(parentEdge)) {
                // edge has already been divided by adjacent face
                childedges = edges.kids(parentEdge);
            }
            else {
                // divide parent edge
                childedges[0] = edges.n();
                childedges[1] = edges.n() + 1;
                edges.divide(parentEdge, particles);
            }
            
            // connect child edge vertex to child faces
            newFaceVerts[i][(i+1)%4] = edges.dest(childedges[0]);
            newFaceVerts[(i+1)%4][i] = edges.dest(childedges[0]);
               
            // connect child edges to child faces
            if (edges.positiveOrientation(parentEdge, ind)) {
                newFaceEdges[i][i] = childedges[0];
                edges.setLeftFace(childedges[0], _faces.size() + i);
                
                newFaceEdges[(i+1)%4][i] = childedges[1];
                edges.setLeftFace(childedges[1], _faces.size() + (i+1)%4);
            }
            else {
                newFaceEdges[i][i] = childedges[1];
                edges.setRightFace(childedges[1], _faces.size() +i);
                
                newFaceEdges[(i+1)%4][i] = childedges[0];
                edges.setRightFace(childedges[0], _faces.size() + (i+1)%4);
            }
        }
        
        std::array<index_type, 4> newFaceInds;
        std::array<index_type, 4> newEdgeInds;
        for (int i=0; i<4; ++i) {
            newFaceInds[i] = _faces.size() + i;
            newEdgeInds[i] = edges.n() + i;
        }
        
        edges.insert(newFaceVerts[0][1], newFaceVerts[0][2], newFaceInds[0], newFaceInds[1]);
        newFaceEdges[0][1] = newEdgeInds[0];
        newFaceEdges[1][3] = newEdgeInds[0];
        
        edges.insert(newFaceVerts[2][0], newFaceVerts[2][3], newFaceInds[3], newFaceInds[2]);
        newFaceEdges[3][1] = newEdgeInds[1];
        newFaceEdges[2][3] = newEdgeInds[1];
        
        edges.insert(newFaceVerts[2][1], newFaceVerts[2][0], newFaceInds[1], newFaceInds[2]); 
        newFaceEdges[2][1] = newEdgeInds[2];
        newFaceEdges[0][2] = newEdgeInds[2];
        
        edges.insert(newFaceVerts[3][1], newFaceVerts[3][0], newFaceInds[0], newFaceInds[3]);
        newFaceEdges[0][3] = newEdgeInds[3];
        newFaceEdges[2][0] = newEdgeInds[3];
        
        // create new center particles
        std::vector<Vec<ndim>> quadCorners(4);
        std::vector<Vec<ndim>> lagQuadCorners(4);
        std::vector<Vec<ndim>> newFaceCenters(4);
        std::vector<Vec<ndim>> newFaceLagCenters(4);
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                quadCorners[j] = particles.physCrd(newFaceVerts[i][j]);
                lagQuadCorners[j] = particles.lagCrd(newFaceVerts[i][j]);
            }   
            if (_geom == SPHERICAL_SURFACE_GEOMETRY) {
                newFaceCenters[i] = sphereBaryCenter(quadCorners, _radius);
                newFaceLagCenters[i] = sphereBaryCenter(lagQuadCorners, _radius);
            }
            else {
                newFaceCenters[i] = baryCenter(quadCorners);
                newFaceLagCenters[i] = baryCenter(lagQuadCorners);
            }
        }
        // create child faces
        const index_type faceInsert = _faces.size();
        for (int i=0; i<4; ++i) {
            particles.insert(newFaceCenters[i], newFaceLagCenters[i]);
            this->insert(std::vector<index_type>(particles.n()-1), newFaceVerts[i], newFaceEdges[i], ind);
            this->_faces[ind]->_kids[i] = faceInsert  + i;
        }
        for (index_type i=faceInsert; _faces.size(); ++i) {
            this->_faces[i]->setArea(_geom, particles, _radius);
        }
    }
    else if (dynamic_cast<QuadCubicFace<ndim>*>(_faces[0].get())) {
        throw std::runtime_error("QuadCubicFaces not implemented.");
    }
    this->_nActive += 3;
    this->_faces[ind]->_area = 0.0;

}


template class FaceSet<2>;
template class FaceSet<3>;
}

