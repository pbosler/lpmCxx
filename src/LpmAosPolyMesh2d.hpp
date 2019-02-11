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

#ifdef HAVE_KOKKOS
#include "Kokkos_Core.hpp"
#include "Kokkos_View.hpp"
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
        
        inline ParticleSet<ndim>* particle_set_raw_ptr() {return _particles.get();}
        
        inline index_type nActive() const {return _particles->nActive();}
        inline index_type nPassive() const {return _particles->nPassive();}
        
        inline Vec<ndim> physCrd(const index_type pind) const {return _particles->physCrd(pind);}
        inline Vec<ndim> lagCrd(const index_type pind) const {return _particles->lagCrd(pind);}
        inline scalar_type weight(const index_type pind) const {return _particles->weight(pind);}
        
        inline bool isPassive(const index_type pind) const {return _particles->isVertex(pind);}
        inline bool isActive(const index_type pind) const {return !_particles->isVertex(pind);}
        inline index_type n() const {return _particles->n();}
        
        void registerScalarField(const std::string& name) {_particles->registerScalarField(name);}
        void registerVectorField(const std::string& name) {_particles->registerVectorField(name);}
        void setScalarFieldValue(const std::string& name, const index_type p_ind, const scalar_type val) {
            _particles->initScalarFieldValue(name, p_ind, val);
        }
        void setVectorFieldValue(const std::string& name, const index_type p_ind, const Vec<ndim>& val) {
            _particles->initVectorFieldValue(name, p_ind, val);
        }
        
#ifdef HAVE_VTK
        vtkSmartPointer<vtkPolyData> toVtkPolyData(const bool useFieldForHeight=false, 
            const std::string field_name=std::string()) const;
#endif
#ifdef HAVE_KOKKOS
        typedef Kokkos::View<scalar_type*[ndim]> vec_view_type;
        typedef typename vec_view_type::HostMirror vec_host_view_type;
        typedef Kokkos::View<scalar_type*> scalar_view_type;
        typedef typename scalar_view_type::HostMirror scalar_host_view_type;
        typedef Kokkos::View<index_type*> index_view_type;
        typedef typename index_view_type::HostMirror index_host_view_type;
        
        /**  Initializes Kokkos::Views for active/passive particles
        */
        void init_pack_coord_views(vec_view_type acrds, vec_host_view_type acrds_h, 
            vec_view_type pcrds, vec_host_view_type pcrds_h, 
            scalar_view_type wgts, scalar_host_view_type wgts_h) const {
            
            acrds = vec_view_type("active_coords", this->nActive());
            pcrds = vec_view_type("passive_coords", this->nPassive());
            wgts = scalar_view_type("weights", this->nActive());
            
            acrds_h = Kokkos::create_mirror_view(acrds);
            pcrds_h = Kokkos::create_mirror_view(pcrds);
            wgts_h = Kokkos::create_mirror_view(wgts);
            // Fill from data on host
            index_type a_ind=0;
            index_type p_ind=0;
            for (index_type i=0; i<_particles->n(); ++i) {
                const Vec<ndim> pc = this->physCrd(i);
                if (this->isActive(i)) {
                    for (short j=0; j<ndim; ++j) {
                        acrds_h(a_ind, j) = pc[j];
                    }
                    wgts_h(a_ind++) = this->weight(i);
                }
                else {
                    for (short j=0; j<ndim; ++j) {
                        pcrds(p_ind,j) = pc[j];
                    }
                    p_ind += 1;
                }
            }
            // Copy data to device
            Kokkos::deep_copy(acrds, acrds_h);
            Kokkos::deep_copy(pcrds, pcrds_h);
            Kokkos::deep_copy(wgts, wgts_h);
        }
        
        void init_pack_active_inds(index_view_type a_inds) const {
            a_inds = index_view_type("active_indices", this->nActive());
            index_host_view_type aih = Kokkos::create_mirror_view(a_inds);
            index_type working_ind = 0;
            for (index_type i=0; i<_particles->n(); ++i) {
                if (this->isActive(i)) {
                    aih[working_ind++] = i;
                }
            }
            Kokkos::deep_copy(a_inds, aih);
        }
        
        void init_pack_all_coords(vec_view_type crds, vec_host_view_type hcrds, 
            scalar_view_type wgts, scalar_host_view_type wgts_h) const {
            
            crds = vec_view_type("phys_coords", _particles->n());
            wgts = scalar_view_type("weights", _particles->n());
            hcrds = Kokkos::create_mirror_view(crds);
            wgts_h = Kokkos::create_mirror_view(wgts);
            for (index_type i=0; i<_particles->n(); ++i) {
                const Vec<ndim> pc = this->physCrd(i);
                for (short j=0; j<ndim; ++j) {
                    hcrds(i,j) = pc[j];
                }
                wgts_h[i] = this->weight(i);
            }
            
            Kokkos::deep_copy(crds, hcrds);
            Kokkos::deep_copy(wgts, wgts_h);
        } 
        
        void init_pack_vector_field(const std::string& name, vec_view_type vf, vec_host_view_type hvf) const {
            vf = vec_view_type(name, _particles->n());
            hvf = Kokkos::create_mirror_view(vf);
            for (index_type i=0; i<_particles->n(); ++i) {
                const Vec<ndim> fval = _particles->vecVal(i, name);
                for (short j=0; j<ndim; ++j) {
                    hvf(i,j) = fval[j];
                }
            }
            Kokkos::deep_copy(vf, hvf);
        }      
        
        void init_pack_scalar_field(const std::string& name, scalar_view_type sf, scalar_host_view_type hsf) const {
            sf = scalar_view_type(name, _particles->n());
            hsf = Kokkos::create_mirror_view(sf);
            for (index_type i=0; i<_particles->n(); ++i) {
                hsf(i) = _particles->scalarVal(i, name);
            }
            Kokkos::deep_copy(sf, hsf);
        }
        
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
