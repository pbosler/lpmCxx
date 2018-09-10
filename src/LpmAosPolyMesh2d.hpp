#ifndef LPM_AOS_POLYMESH_2D_HPP
#define LPM_AOS_POLYMESH_2D_HPP

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosTypes.hpp"
#include "LpmAosParticle.hpp"
#include "LpmAosParticleFactory.hpp"
#include "LpmAosParticleSet.hpp"
#include "LpmAosEdge.hpp"
#include "LpmAosEdgeFactory.hpp"
#include "LpmAosEdgeSet.hpp"
#include "LpmAosFace.hpp"
#include "LpmAosFaceFactory.hpp"
#include "LpmAosFaceSet.hpp"
#include "LpmAosMeshSeed.hpp"

#ifdef HAVE_VTK
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyData.h"
#endif

namespace Lpm {
namespace Aos {

///  2D meshes made of vertices, edges, faces; possibly embedded in R3.
/*
    All data are carried by particles; there are particles associated with vertices always, and there may be 
    other particles associated with edges or faces (for variable staggering, e.g.).
    
    EdgeSet and FaceSet record connectivity only.
*/
template <int ndim> class PolyMesh2d {
    public:
        PolyMesh2d(const std::shared_ptr<MeshSeed<ndim>> seed, const std::shared_ptr<ParticleFactory<ndim>> pfac, 
            const std::shared_ptr<EdgeFactory<ndim>> efac, const std::shared_ptr<FaceFactory<ndim>> ffac, 
            const index_type initnest, const index_type maxnest, const index_type amrlimit, 
            const scalar_type domainRadius=1.0, const scalar_type t=0.0, const index_type tind=0) : 
            _seed(seed), _pfac(pfac), _efac(efac), _ffac(ffac), _initnest(initnest), _maxnest(maxnest), _amrlimit(amrlimit),
            _radius(domainRadius), _time(t), _time_index(tind)
            {}
        
        void initStaggeredVerticesAndFacesFromSeed();

        std::string infoString(const bool printAll=false) const;
        
        inline index_type nLeafFaces() const {return _faces->nLeaves();}
        inline index_type nLeafEdges() const {return _edges->nLeaves();}
        inline index_type nVertices() const {return _particles->n() - _faces->nIntrsPerFace()*_faces->nLeaves();}
        
        inline scalar_type surfaceArea() const {return _faces->totalArea();}
        
#ifdef HAVE_VTK
        vtkSmartPointer<vtkPolyData> toVtkPolyData(const bool useFieldForHeight=false, 
            const std::string field_name=std::string()) const;
#endif
    
    protected:
        short _initnest;
        short _maxnest;
        short _amrlimit;
        scalar_type _radius;
        
        scalar_type _time;
        index_type _time_index;
        
#ifdef HAVE_VTK        
//         vtkSmartPointer<vtkPoints> verticesToVtkPoints(const bool useFieldForHeight=false, 
//             const std::string scalarFieldName=std::string()) const;
        vtkSmartPointer<vtkPointData> fieldsToVtkPointData() const;
        vtkSmartPointer<vtkCellData> fieldsToVtkCellData() const;
#endif
        
        std::shared_ptr<MeshSeed<ndim>> _seed;
        
        std::shared_ptr<ParticleFactory<ndim>> _pfac;
        std::shared_ptr<EdgeFactory<ndim>> _efac;
        std::shared_ptr<FaceFactory<ndim>> _ffac;
        
        std::unique_ptr<ParticleSet<ndim>> _particles;
        std::unique_ptr<EdgeSet<ndim>> _edges;
        std::unique_ptr<FaceSet<ndim>> _faces;   
};

}
}
#endif
