#include "LpmAosPolyMesh2d.hpp"
#ifdef HAVE_VTK
#include "vtkDoubleArray.h"
#endif

namespace Lpm {
namespace Aos {

template <int ndim> PolyMesh2d<ndim>::PolyMesh2d<ndim>(const MeshSeed* seed, const index_type initRecursion, 
    const index_type maxRecursion, const int amrLimit, const scalar_type domainRadius=1.0) :
    seedId(seed->idString()), _initRecursion(initRecursion), _maxRecursion(maxRecursion), _amrLimit(amrLimit), 
    _radius(domainRadius), _tindex(0), _time(0.0) {
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
            vals->SetName(sfields[i]);
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
            vals->SetName(vfields[i]);
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
