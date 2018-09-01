#include "LpmAosPolyMesh2d.hpp"
#ifdef HAVE_VTK
#include "vtkDoubleArray.h"
#endif

namespace Lpm {
namespace Aos {

template <int ndim> void PolyMesh2d<ndim>::initFromSeedStaggeredFacesAndVerts(const MeshSeed* seedptr, 
    const std::shared_ptr<ParticleFactory<ndim>> pfac, const std::shared_ptr<EdgeFactory<ndim>> efac,
    const std::shared_ptr<FaceFactory<ndim>> ffac) {
    _seedId = seedptr->idString();
    _geom = seedptr->geometryType();
    
    index_type nMaxVerts;
    index_type nMaxFaces;
    index_type nMaxEdges;
    seedptr->determineMaxAllocations(nMaxVerts, nMaxEdges, nMaxFaces, _maxRecursion);
    
    _vertexParticles = std::unique_ptr<ParticleSet<ndim>>(new ParticleSet<ndim>(pfac, nMaxVerts));
    _faceParticles = std::unique_ptr<ParticleSet<ndim>>(new ParticleSet<ndim>(pfac, nMaxFaces));
    
    _edges = std::unique_ptr<EdgeSet<ndim>>(new EdgeSet<ndim>(efac, _geom, nMaxEdges));
    _faces = std::unique_ptr<FaceSet<ndim>>(new FaceSet<ndim>(ffac, _geom, nMaxFaces, _radius));
    
    for (index_type i=0; i<seedptr->nCrds(); ++i) {
        switch (ndim) {
            case (2) : {
                _vertexParticles->insert(seedptr->r2Crds[i]);
                break;
            }
            case (3) : {
                _vertexParticles->insert(seedptr->r3Crds[i]);
                break;
            }
        }
    }
    if (seedptr->faceCrdsIncluded()) {
        for (index_type i=0; i<seedptr->nFaces(); ++i) {
            switch (ndim) {
                case (2) : {
                    _faceParticles->insert(seedptr->r2Crds[seedptr->nVerts()+i]);
                    break;
                }
                case (3) : {
                    _faceParticles->insert(seedptr->r3Crds[seedptr->nVerts()+i]);
                    break;
                }
            }
        }
    }
    
    for (index_type i=0; i<seedptr->nEdges(); ++i) {
        _edges->insert(seedptr->edgeOrigs[i], seedptr->edgeDests[i], seedptr->edgeLefts[i], seedptr->edgeRights[i]);
        _edges->enrich(i, *_vertexParticles);
    }
    const std::vector<index_type> intrs(1,-1);
    const index_type root_face_parent_index = -1;
    for (index_type i=0; i<seedptr->nFaces(); ++i) {
        intrs[0] = i;
        _faces->insert(intrs, seedptr->faceVerts[i], seedptr->faceEdges[i], root_face_parent_index);
        _faces->enrich(i, *_vertexParticles, *_edges);
    }
}



#ifdef HAVE_VTK
template <int ndim>  vtkSmartPointer<vtkPolyData> PolyMesh2d<ndim>::toVtkPolyData() const {
        vtkSmartPointer<vtkPolyData> result = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> pts = _vertexParticles->toVtkPoints();
        vtkSmartPointer<vtkCellArray> polys = _faces->toVtkCellArray();
        result->SetPoints(pts);
        result->SetPolys(polys);
        const std::vector<std::string> sfields = _vertexParticles->getScalarFieldNames();
        for (int i=0; i<sfields.size(); ++i) {
            vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
            vals->SetName(sfields[i].c_str());
            vals->SetNumberOfComponents(1);
            vals->SetNumberOfTuples(_vertexParticles->n());
            for (int j=0; j<_vertexParticles.n(); ++j) {
                vals->InsertTuple1(_vertexParticles->scalarVal(j, sfields[i]));
            }
            result->GetPointData()->AddArray(vals);
        }
        const std::vector<std::string> vfields = _vertexParticles->getVectorFieldNames();
        for (int i=0; i<vfields.size(); ++i) {
            vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();            
            vals->SetName(vfields[i].c_str());
            vals->SetNumberOfComponents(ndim);
            vals->SetNumberOfTuples(_vertexParticles.n());
            switch (ndim) {
                case (2) : {
                    for (index_type j=0; j<_vertexParticles.n(); ++j) {
                        vals->InsertTuple2(_vertexParticles->vectorVal(j, vfields[i]));
                    }
                    break;
                }
                case (3) : {
                    for (index_type j=0; j<_vertexParticles.n(); ++j) {
                        vals->InsertTuple3(_vertexParticles->vectorVal(j, vfields[i]));
                    }
                    break;
                }
            }
            result->GetPointData()->AddArray(vals);
        }
        return result;
    }
#endif

}
}
