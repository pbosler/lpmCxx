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
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyData.h"
#endif

namespace Lpm {
namespace Aos {

template <int ndim> class PolyMesh2d {
    public:
        PolyMesh2d(const MeshSeed* seed, const std::shared_ptr<ParticleFactory<ndim>> pfac,
            const std::shared_ptr<EdgeFactory<ndim>> efac, const std::shared_ptr<FaceFactory<ndim>> ffac, 
            const index_type nMaxParticles, const index_type nMaxEdges, const index_type nMaxFaces,
            const index_type initRecursion, const index_type maxRecursion, 
            const int amrLimit, const scalar_type domainRadius=1.0) : 
            _seed(seed), _pfac(pfac), _efac(efac), _ffac(ffac), 
            _initRecursion(initRecursion), _maxRecursion(maxRecursion), _amrLimit(amrLimit), 
            _radius(domainRadius), _tindex(0), _time(0.0), 
            _vertexParticles(pfac, nMaxParticles),
            _edges(efac, seed->geometryType(), nMaxEdges),
            _faces(ffac, nMaxFaces, seed->geometryType(), domainRadius)
            {}
        
 //        void initFromSeedStaggeredFacesAndVerts(const MeshSeed* seedptr);
//         
// //         todo
//         void initFromSeedCollocated(const MeshSeed* seedptr, const std::shared_ptr<ParticleFactory<ndim>> pfac,
//             const std::shared_ptr<EdgeFactory<ndim>> efac, const std::shared_ptr<FaceFactory<ndim>> ffac) {};
//         
// //          todo
//         void initFromFileCollocated(const std::string filename){};
//         
// //         todo
//         void initFromSeedStaggeredFacesEdgesAndVerts(const MeshSeed* seedptr, const std::shared_ptr<ParticleFactory<ndim>> pfac,
//             const std::shared_ptr<EdgeFactory<ndim>> efac, const std::shared_ptr<FaceFactory<ndim>> ffac) {};
        
        inline GeometryType geometryType() const {return _geom;}
        
//         ParticleSet<ndim>* verticesPtr() {return _vertexParticles.get();}
//         ParticleSet<ndim>* edgeParticlePtr() {return _edgeParticles.get();}
//         ParticleSet<ndim>* faceParticlePtr() {return _faceParticles.get();}
//         EdgeSet<ndim>* edgesPtr() {return _edges.get();}
//         FaceSet<ndim>* facesPtr() {return _faces.get();}

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
        GeometryType _geom;
        FaceType _facekind;
        const MeshSeed* _seed;
        
        bool verifyFactories(const MeshSeed* seed, const std::shared_ptr<FaceFactory<ndim>> ffac);
        
        void determineMaxAllocations(index_type& nv, index_type& nf, index_type& ne) const;
    
        std::shared_ptr<ParticleFactory<ndim>> _pfac;
        std::shared_ptr<EdgeFactory<ndim>> _efac;
        std::shared_ptr<FaceFactory<ndim>> _ffac;
    
        ParticleSet<ndim> _vertexParticles;
        ParticleSet<ndim> _edgeParticles;
        ParticleSet<ndim> _faceParticles;
        EdgeSet<ndim> _edges;
        FaceSet<ndim> _faces;
        
        index_type nRootFaces;
        
//         index_type walkSearch(const Vec<ndim>& querypt, const index_type startIndex) const;
//         index_type treeSearch(const Vec<ndim>& querypt, const index_type, startIndex) const;
//         index_type nearestRootFace(const Vec<ndim>& querypt) const;
    
};

}
}
#endif
