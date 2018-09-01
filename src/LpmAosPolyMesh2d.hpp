#ifndef LPM_AOS_POLYMESH_2D_HPP
#define LPM_AOS_POLYMESH_2D_HPP

#include <vector>
#include <memory>
#include <string>
#include <exception>
#include <algorithm>
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
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyData.h"
#endif

namespace Lpm {
namespace Aos {

template <int ndim> class PolyMesh2d {
    public:
        PolyMesh2d(const MeshSeed* seed, const index_type initRecursion, const index_type maxRecursion, 
            const int amrLimit, const scalar_type domainRadius=1.0);
        
        inline GeometryType geom() const {return _faces->geom;}
        
        ParticleSet<ndim>* verticesPtr() {return _vertexParticles.get();}
        ParticleSet<ndim>* edgeParticlePtr() {return _edgeParticles.get();}
        ParticleSet<ndim>* faceParticlePtr() {return _faceParticles.get();}
        EdgeSet<ndim>* edgesPtr() {return _edges.get();}
        FaceSet<ndim>* facesPtr() {return _faces.get();}

#ifdef HAVE_VTK
        vtkSmartPointer<vtkPolyData> toVtkPolyData() const;
#endif
    
        void writePolydataFile(std::ostream& os) const;        
        
    protected:
        std::string _seedId;
        int _amrLimit;
        int _initRecursion;
        int _maxRecursion;
        scalar_type _radius;
        index_type _tindex;
        scalar_type _time;
    
        std::unique_ptr<ParticleSet<ndim>> _vertexParticles;
        std::unique_ptr<ParticleSet<ndim>> _edgeParticles;
        std::unique_ptr<ParticleSet<ndim>> _faceParticles;
        std::unique_ptr<EdgeSet<ndim>> _edges;
        std::unique_ptr<FaceSet<ndim>> _faces;
        
        index_type nRootFaces;
        
        index_type walkSearch(const Vec<ndim>& querypt, const index_type startIndex) const;
        index_type treeSearch(const Vec<ndim>& querypt, const index_type, startIndex) const;
        index_type nearestRootFace(const Vec<ndim>& querypt) const;
        
        void initFromSeed(const MeshSeed& seed, const index_type maxRecursion, const scalar_type radius);
    
};

}
}
#endif
