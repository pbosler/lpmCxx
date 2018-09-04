#include "LpmAosPolyMesh2d.hpp"
#include <exception>

#ifdef HAVE_VTK
#include "vtkDoubleArray.h"
#endif

namespace Lpm {
namespace Aos {

template <int ndim> void PolyMesh2d<ndim>::initFromSeedStaggeredFacesAndVerts() {
    std::shared_ptr<MeshSeed<ndim>> seed = _seed.lock();
    const index_type nseedverts = seed->nVerts();
    const index_type nseededges = seed->nEdges();
    const index_type nseedfaces = seed->nFaces();
    const index_type nseedcrds = seed->nCrds();
    const int nedgesperface = seed->nEdgesPerFace();
    
    index_type nMaxVertices;
    index_type nMaxEdges;
    index_type nMaxFaces;
    seed->determineMaxAllocations(nMaxVertices, nMaxEdges, nMaxFaces, _maxRecursion);
    
    std::cout << "PolyMesh2d::initFromSeedStaggeredFacesAndVerts using meshSeed " << seed->idString() << "; reserving space for "
    << nMaxVertices << " vertices, " << nMaxEdges << " edges, and " << nMaxFaces << " faces." << std::endl;
    
    _vertexParticles = std::unique_ptr<ParticleSet<ndim>>(new ParticleSet<ndim>(_pfac, nMaxVertices));
    _faceParticles = std::unique_ptr<ParticleSet<ndim>>(new ParticleSet<ndim>(_pfac, nMaxFaces));
    _edges = std::unique_ptr<EdgeSet<ndim>>(new EdgeSet<ndim>(_efac, _geom, nMaxEdges));
    _faces = std::unique_ptr<FaceSet<ndim>>(new FaceSet<ndim>(_ffac, nMaxFaces, _geom, _radius));
    
//     std::cout << _edges->infoString();

    for (index_type i=0; i<nseedverts; ++i) {
        _vertexParticles->insert(seed->crds[i]);
    }
    
    if (seed->faceCrdsIncluded()) {
        for (index_type i=0; i<nseedfaces; ++i) {
            _faceParticles->insert(seed->crds[i+nseedverts]);    
        }
    }
    else {
        std::vector<Vec<ndim>> verts(nedgesperface);
        for (index_type i=0; i<nseedfaces; ++i) {
            for (index_type j=0; j<nedgesperface; ++j) {
                verts[j] = seed->crds[seed->faceVerts[i][j]];
            }
            if (_geom == SPHERICAL_SURFACE_GEOMETRY ){
                _faceParticles->insert(sphereBaryCenter(verts));
            }
            else {
                _faceParticles->insert(baryCenter(verts));
            }
        }
    }

    
    for (index_type i=0; i<nseededges; ++i) {
        _edges->insert(seed->edgeOrigs[i], seed->edgeDests[i], seed->edgeLefts[i], seed->edgeRights[i]);
    }
    std::vector<index_type> interior_ind(1,-1);
    const index_type root_parent_ind = -1;
    for (index_type i=0; i<nseedfaces; ++i) {
        interior_ind[0] = i;
        _faces->insert(interior_ind, seed->faceVerts[i], seed->faceEdges[i], root_parent_ind, _radius);
    }

    std::cout << "PolyMesh2d::initFromSeedStaggeredFacesAndVerts: starting uniform refinement." << std::endl;
}


// template <int ndim> void PolyMesh2d<ndim>::initFromSeedStaggeredFacesAndVerts(const MeshSeed* seedptr, 
//      const std::shared_ptr<ParticleFactory<ndim>> pfac,  const std::shared_ptr<EdgeFactory<ndim>> efac,
//      const std::shared_ptr<FaceFactory<ndim>> ffac) 
// {
//     _seedId = seedptr->idString();
//     _geom = seedptr->geometryType();
//     if (!verifyFactories(seedptr, ffac)) {
//         throw std::runtime_error("mismatched face type.");
//     }
//     _facekind = seedptr->faceType();
//     
//     index_type nMaxVerts;
//     index_type nMaxFaces;
//     index_type nMaxEdges;
//     seedptr->determineMaxAllocations(nMaxVerts, nMaxEdges, nMaxFaces, _maxRecursion);
//     
//     _vertexParticles = ParticleSet<ndim>(pfac, nMaxVerts);
//     _faceParticles = ParticleSet<ndim>(pfac, nMaxFaces);
//     //_edgeParticles = nullptr;
//     
//     _edges = EdgeSet<ndim>(efac, _geom, nMaxEdges);
//     _faces = FaceSet<ndim>(ffac, nMaxFaces, _geom, _radius);
//     
//     for (index_type i=0; i<seedptr->nCrds(); ++i) {
//         switch (ndim) {
//             case (2) : {
//                 _vertexParticles.insert(seedptr->r2Crds[i]);
//                 break;
//             }
//             case (3) : {
//                 _vertexParticles.insert(seedptr->r3Crds[i]);
//                 break;
//             }
//         }
//     }
//     if (seedptr->faceCrdsIncluded()) {
//         for (index_type i=0; i<seedptr->nFaces(); ++i) {
//             switch (ndim) {
//                 case (2) : {
//                     _faceParticles.insert(seedptr->r2Crds[seedptr->nVerts()+i]);
//                     break;
//                 }
//                 case (3) : {
//                     _faceParticles.insert(seedptr->r3Crds[seedptr->nVerts()+i]);
//                     break;
//                 }
//             }
//         }
//     }
//     else {
//         for (index_type i=0; i<seedptr->nFaces(); ++i) {
//             switch (ndim) {
//                 case (2) : {
//                     std::vector<Vec<2>> verts(seedptr->nEdgesPerFace());
//                     for (int j=0; j<seedptr->nEdgesPerFace(); ++j) {
//                         verts[j] = seedptr->r2Crds[seedptr->faceVerts[i][j]];    
//                     }
//                     _faceParticles.insert(baryCenter(verts));
//                     break;
//                 }
//                 case (3) : {
//                     std::vector<Vec<3>> verts(seedptr->nEdgesPerFace());
//                     for (int j=0; j<seedptr->nEdgesPerFace(); ++j) {
//                         verts[j] = seedptr->r3Crds[seedptr->faceVerts[i][j]];
//                     }
//                     if (_geom == SPHERICAL_SURFACE_GEOMETRY) {
//                         _faceParticles.insert(sphereBaryCenter(verts));
//                     }
//                     else {
//                         _faceParticles.insert(baryCenter(verts));
//                     }
//                     break;
//                 }
//             }
//         }
//     }
//     
//     for (index_type i=0; i<seedptr->nEdges(); ++i) {
//         _edges.insert(seedptr->edgeOrigs[i], seedptr->edgeDests[i], seedptr->edgeLefts[i], seedptr->edgeRights[i]);
//         _edges.enrich(i, _vertexParticles);
//     }
//     std::vector<index_type> intrs(1,-1);
//     const index_type root_face_parent_index = -1;
//     for (index_type i=0; i<seedptr->nFaces(); ++i) {
//         intrs[0] = i;
//         _faces.insert(intrs, seedptr->faceVerts[i], seedptr->faceEdges[i], root_face_parent_index);
//         _faces.enrich(i, _vertexParticles, _edges);
//     }
// }

template <int ndim> bool PolyMesh2d<ndim>::verifyFactories(const MeshSeed<ndim>* seed, const std::shared_ptr<FaceFactory<ndim>> ffac) {    
    return (seed->faceType() == ffac->basicFaceType());
}


#ifdef HAVE_VTK
template <int ndim>  vtkSmartPointer<vtkPolyData> PolyMesh2d<ndim>::toVtkPolyData() const {
        vtkSmartPointer<vtkPolyData> result = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> pts = _vertexParticles->toVtkPoints();
        vtkSmartPointer<vtkCellArray> polys = _faces->toVtkCellArray();
        vtkSmartPointer<vtkPointData> ptdata = _vertexParticles->fieldsToVtkPointData();
        vtkSmartPointer<vtkCellData> celldata = _faceParticles->fieldsToVtkCellData();
        result->SetPoints(pts);
        result->SetPolys(polys);
        const int nPtFields = ptdata->GetNumberOfArrays();
        const int nCellFields = celldata->GetNumberOfArrays();
        for (int i=0; i<nPtFields; ++i) {
            result->GetPointData()->AddArray(ptdata->GetAbstractArray(i));
        }
        for (int i=0; i<nCellFields; ++i) {
            result->GetCellData()->AddArray(celldata->GetAbstractArray(i));
        }
//         std::cout << "\tmesh data written." << std::endl;
//         const std::vector<std::string> sfields = _vertexParticles->getScalarFieldNames();
//         std::cout << "\tfound " << sfields.size() << " scalar fields." << std::endl;
//         for (int i=0; i<sfields.size(); ++i) std::cout << sfields[i] << std::endl;
//         for (int i=0; i<sfields.size(); ++i) {
//             std::cout << "\twriting field: " << sfields[i] << std::endl;
//             vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
//             vals->SetName(sfields[i].c_str());
//             vals->SetNumberOfComponents(1);
//             vals->SetNumberOfTuples(_vertexParticles->n());
//             for (index_type j=0; j<_vertexParticles->n(); ++j) {
//                 vals->InsertTuple1(j, _vertexParticles->scalarVal(j, sfields[i]));
//             }
//             result->GetPointData()->AddArray(vals);
//         }
//         const std::vector<std::string> vfields = _vertexParticles->getVectorFieldNames();
//         for (int i=0; i<vfields.size(); ++i) {
//             vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();            
//             vals->SetName(vfields[i].c_str());
//             vals->SetNumberOfComponents(ndim);
//             vals->SetNumberOfTuples(_vertexParticles->n());
//             std::cout << "\twriting field: " << vfields[i] << std::endl;
//             switch (ndim) {
//                 case (2) : {
//                     for (index_type j=0; j<_vertexParticles->n(); ++j) {
//                         const std::vector<scalar_type> vecval = _vertexParticles->vectorVal(j, vfields[i]);
//                         vals->InsertTuple2(j, vecval[0], vecval[1]);
//                     }
//                     break;
//                 }
//                 case (3) : {
//                     for (index_type j=0; j<_vertexParticles->n(); ++j) {
//                         const std::vector<scalar_type> vecval = _vertexParticles->vectorVal(j, vfields[i]);
//                         vals->InsertTuple3(j, vecval[0], vecval[1], vecval[2]);
//                     }
//                     break;
//                 }
//             }
//             result->GetPointData()->AddArray(vals);
//         }
        return result;
    }
#endif


template class PolyMesh2d<2>;
template class PolyMesh2d<3>;
}
}
