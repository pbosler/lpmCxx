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
            const index_type initRecursion, const index_type maxRecursion, 
            const int amrLimit, const scalar_type domainRadius=1.0) : 
            _seed(seed), _pfac(pfac), _efac(efac), _ffac(ffac), 
            _initRecursion(initRecursion), _maxRecursion(maxRecursion), _amrLimit(amrLimit), 
            _radius(domainRadius), _tindex(0), _time(0.0), 
            _seedId(seed->idString()), _geom(seed->geometryType()), _facekind(seed->faceType())
            {}
        
        void initFromSeedStaggeredFacesAndVerts();
        
        inline index_type nParticles() const {return _vertexParticles->n();}
        inline index_type nMaxParticles() const {return _vertexParticles->nMax();}
        
        inline index_type nEdges() const {return _edges->n();}
        inline index_type nLeafEdges() const {return _edges->nLeaves();}
        inline index_type nMaxEdges() const {return _edges->nMax();}
        
        inline index_type nFaces() const {return _faces->n();}
        inline index_type nLeafFaces() const {return _faces->nLeaves();}
        inline index_type nMaxFaces() const {return _faces->nMax();}
        
        
        
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
        
        std::string infoString() const;

#ifdef HAVE_VTK
        vtkSmartPointer<vtkPolyData> toVtkPolyData() const;
#endif

    // ascii
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
        std::weak_ptr<MeshSeed<ndim>> _seed;
        
        bool verifyFactories(const MeshSeed<ndim>* seed, const std::shared_ptr<FaceFactory<ndim>> ffac);
        
        void determineMaxAllocations(index_type& nv, index_type& nf, index_type& ne) const;
    
        std::shared_ptr<ParticleFactory<ndim>> _pfac;
        std::shared_ptr<ParticleFactory<ndim>> _epfac;
        std::shared_ptr<ParticleFactorY<ndim>> _fpfac;
        std::shared_ptr<EdgeFactory<ndim>> _efac;
        std::shared_ptr<FaceFactory<ndim>> _ffac;
    
        std::unique_ptr<ParticleSet<ndim>> _vertices;
        std::unique_ptr<EdgeSet<ndim>> _edges;
        std::unique_ptr<FaceSet<ndim>> _faces;
        
        index_type nRootFaces;
        
//         index_type walkSearch(const Vec<ndim>& querypt, const index_type startIndex) const;
//         index_type treeSearch(const Vec<ndim>& querypt, const index_type, startIndex) const;
//         index_type nearestRootFace(const Vec<ndim>& querypt) const;
    
};

}
}
#endif
