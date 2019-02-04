#include "LpmAosPolyMesh2d.hpp"
#include <exception>

#ifdef HAVE_VTK
#include "vtkDoubleArray.h"
#endif

namespace Lpm {
namespace Aos {

template <int ndim> void PolyMesh2d<ndim>::initStaggeredVerticesAndFacesFromSeed() {
    // determine memory allocations
    index_type nmaxparticles;
    index_type nmaxedges;
    index_type nmaxfaces;
    _seed->determineMaxAllocations(nmaxparticles, nmaxedges, nmaxfaces, _maxnest);
    std::cout << "polymesh memory requirements: " << nmaxparticles << " vertices, " << nmaxedges << " edges, " << nmaxfaces << " faces." << std::endl;
    // build basic containers
    const GeometryType geom = _seed->geometryType();
    _particles = std::unique_ptr<ParticleSet<ndim>>(new ParticleSet<ndim>(_pfac, nmaxparticles+nmaxfaces));
    _edges = std::unique_ptr<EdgeSet<ndim>>(new EdgeSet<ndim>(_efac, geom, nmaxedges));
    _faces = std::unique_ptr<FaceSet<ndim>>(new FaceSet<ndim>(_ffac, nmaxfaces, geom, _radius));
    // build root mesh from seed
    const int nseedverts = _seed->nVerts();
    const int nseedfaces = _seed->nFaces();
    const int nseededges = _seed->nEdges();
    // insert particles from seed coordinates
    for (int i=0; i<nseedverts; ++i) {
        _particles->insert(_seed->crds[i], 0.0, true);
    }
    if (_seed->faceCrdsIncluded()) {
        for (int i=0; i<nseedfaces; ++i) {
            _particles->insert(_seed->crds[nseedverts + i]);
        }
    }
    else {
        const short nvertsperface = _seed->faceVerts[0].size();
        std::vector<Vec<ndim>> corners(nvertsperface);
        for (int i=0; i<nseedfaces; ++i){
            for (short j=0; j<nvertsperface; ++j) {
                corners[j] = _seed->crds[_seed->faceVerts[i][j]];
            }
            if (geom == PLANAR_GEOMETRY || geom == CARTESIAN_3D_GEOMETRY) {
                _particles->insert(baryCenter(corners));
            }
            else if (geom == SPHERICAL_SURFACE_GEOMETRY) {
                _particles->insert(sphereBaryCenter(corners, _radius));
            }
        }
    }
    // insert edges from seed
    for (int i=0; i<nseededges; ++i) {
        _edges->insert(_seed->edgeOrigs[i], _seed->edgeDests[i], _seed->edgeLefts[i], _seed->edgeRights[i],
             _seed->edgeInteriors[i]);
    }
    // insert faces from seed
    const index_type root_parent = -1;
    for (int i=0; i<nseedfaces; ++i) {
        _faces->insert(_seed->faceInteriors[i], _seed->faceVerts[i], _seed->faceEdges[i], root_parent);
    }
    _faces->setArea(*_particles);
    
    // refine to desired uniform level
    index_type start_index = 0;
    for (index_type i=0; i<_initnest; ++i) {
        const index_type nfaces_old = _faces->n();
        for (index_type j=start_index; j<nfaces_old; ++j) {
            _faces->divide(j, *_particles, *_edges);
        }
        start_index = nfaces_old;
    }
}

template <int ndim> std::string PolyMesh2d<ndim>::infoString(const bool printAll) const {
    std::ostringstream ss;
    ss << "POLYMESH INFO" << std::endl;
    ss << "seed id: " << _seed->idString() << std::endl;
    ss << "initnest = " << _initnest << std::endl;
    ss << "maxnest = " << _maxnest << std::endl;
    ss << "amrlimit = " << _amrlimit << std::endl;
ss << "radius = " << _radius << std::endl;
    ss << _particles->infoString(printAll);
    ss << _edges->infoString(printAll);
    ss << _faces->infoString(printAll);
    return ss.str();
}

#ifdef HAVE_VTK
template <int ndim> vtkSmartPointer<vtkPointData> PolyMesh2d<ndim>::fieldsToVtkPointData() const {
    vtkSmartPointer<vtkPointData> result = vtkSmartPointer<vtkPointData>::New();
//     add geometric weight data
    const index_type nverts = this->nVertices();
    index_type ctr=0;
    // collect field names
    const std::vector<std::string> sfields = _particles->getScalarFieldNames();
    const std::vector<std::string> vfields = _particles->getVectorFieldNames();
    // add field data
    for (int i=0; i<sfields.size(); ++i) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName(sfields[i].c_str());
        data->SetNumberOfComponents(1);
        data->SetNumberOfTuples(nverts);
        ctr=0;
        for (index_type j=0; j<_particles->n(); ++j) {
            if (_particles->isVertex(j)) {
                data->InsertTuple1(ctr++, _particles->scalarVal(j, sfields[i]));
            }
        }
        result->AddArray(data);
    }
    for (int i=0; i<vfields.size(); ++i) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName(vfields[i].c_str());
        data->SetNumberOfComponents(ndim);
        data->SetNumberOfTuples(nVertices());
        ctr=0;
        for (index_type j=0; j<_particles->n(); ++j) {
            if (_particles->isVertex(j)) {
                const std::vector<scalar_type> val = _particles->vectorVal(j, vfields[i]);
                   if (ndim ==2 ) {
                        data->InsertTuple2(ctr++, val[0], val[1]);
                   }
                   else {
                        data->InsertTuple3(ctr++, val[0], val[1], val[2]);
                   }
            }
        }
    }
    return result;
}

template <int ndim> vtkSmartPointer<vtkCellData> PolyMesh2d<ndim>::fieldsToVtkCellData() const {
    vtkSmartPointer<vtkCellData> result = vtkSmartPointer<vtkCellData>::New();
    // add geometric data
    vtkSmartPointer<vtkDoubleArray> area = vtkSmartPointer<vtkDoubleArray>::New();
    area->SetName("area");
    area->SetNumberOfComponents(1);
    area->SetNumberOfTuples(_faces->nActive());
    index_type ctr=0;
    for (index_type i=0; i<_faces->n(); ++i) {
        if (_faces->isLeaf(i)) {
            area->InsertTuple1(ctr++, _faces->area(i));
        }
    }
    result->AddArray(area);
    // collect field names
    const std::vector<std::string> sfields = _particles->getScalarFieldNames();
    const std::vector<std::string> vfields = _particles->getVectorFieldNames();
    const index_type nintrs = _faces->nIntrsPerFace();
    for (int i=0; i<sfields.size(); ++i) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName(sfields[i].c_str());
        data->SetNumberOfComponents(1);
        data->SetNumberOfTuples(_faces->nActive());
        ctr = 0;
        for (index_type j=0; j<_faces->n(); ++j) {
            if (_faces->isLeaf(j)) {
                scalar_type avg = 0.0;
                const std::vector<index_type> int_inds = _faces->interiors(j);    
                for (index_type k=0; k<int_inds.size(); ++k) {
                    avg += _particles->scalarVal(int_inds[k], sfields[i]);
                }
                avg /= nintrs;
                data->InsertTuple1(ctr++, avg);
            }
        }
        result->AddArray(data);
    }
    for (int i=0; i<vfields.size(); ++i) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName(vfields[i].c_str());
        data->SetNumberOfComponents(ndim);
        data->SetNumberOfTuples(_faces->nActive());
        ctr = 0;
        for (index_type j=0; j<_faces->n(); ++j) {
            if (_faces->isLeaf(j)) {
                Vec<ndim> avg;
                const std::vector<index_type> int_inds = _faces->interiors(j);    
                for (index_type k=0; k<int_inds.size(); ++k) {
                    avg += Vec<ndim>(_particles->vectorVal(int_inds[k], vfields[i]));
                }
                avg.scaleInPlace(1.0/nintrs);
                if (ndim == 2) {
                    data->InsertTuple2(ctr++, avg.x[0], avg.x[1]);
                }
                else if (ndim == 3) {
                    data->InsertTuple3(ctr++, avg.x[0], avg.x[1], avg.x[2]);
                }
            }
        }
        result->AddArray(data);
    }
    return result;
}

template <int ndim>  vtkSmartPointer<vtkPolyData> PolyMesh2d<ndim>::toVtkPolyData(const bool useFieldForHeight,
    const std::string field_name) const {
        vtkSmartPointer<vtkPolyData> result = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
        std::vector<index_type> poffset(_particles->n(), 0);
        index_type vert_ctr=0;
        index_type pt_ctr=0;
        for (index_type i=0; i<_particles->n(); ++i) {
            if (_particles->isVertex(i)) {
                const Vec<ndim> pos = _particles->physCrd(i);
                if (ndim==2) {
                    pts->InsertNextPoint(pos.x[0], pos.x[1], 
                        (useFieldForHeight ? _particles->scalarVal(i, field_name) : 0.0));
                }
                else {
                    pts->InsertNextPoint(pos.x[0], pos.x[1], pos.x[2]);
                }
                poffset[i] = pt_ctr;
            }
            else {
                poffset[i] = ++pt_ctr;

            }
        }
        vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
        const index_type ptsPerFace = _faces->nVertsPerFace();
        for (index_type i=0; i<_faces->n(); ++i) {
            if (_faces->isLeaf(i)) {
                std::vector<index_type> verts = _faces->vertices(i);
                for (int j=0; j<ptsPerFace; ++j) {
                    verts[j] -= poffset[verts[j]];
                }
                polys->InsertNextCell(ptsPerFace);
                for (int j=0; j<ptsPerFace; ++j) {
                    polys->InsertCellPoint(verts[j]);
                }
            }
        }
        vtkSmartPointer<vtkPointData> ptdata = this->fieldsToVtkPointData();
        vtkSmartPointer<vtkCellData> celldata = this->fieldsToVtkCellData();
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
        return result;
    }
#endif


template class PolyMesh2d<2>;
template class PolyMesh2d<3>;
}
}
