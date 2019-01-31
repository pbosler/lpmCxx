#ifndef _LPM_AOS_PARTICLE_SET_HPP
#define _LPM_AOS_PARTICLE_SET_HPP

#include <vector>
#include <memory>
#include <string>
#include <exception>
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmAosParticleFactory.hpp"
#include "LpmAnalyticFunctions.h"

#ifdef HAVE_VTK
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#endif

#ifdef HAVE_KOKKOS
#include "Kokkos_Core.hpp"
#include "Kokkos_View.hpp"
#endif

namespace Lpm {
namespace Aos {

/// A ParticleSet is a collection of Particles arranged linearly in memory.
template <int ndim> class ParticleSet {
    public:
    	/// Constructor. Initializes an empty set.
        ParticleSet() : _factory(nullptr), _nMax(0), _nActive(0) {}
    	/// Constructor. Initializes an empty set.
        ParticleSet(const std::shared_ptr<ParticleFactory<ndim>> factory, const index_type nMax) : 
            _factory(factory), _nMax(nMax), _nActive(0) {
            _particles.reserve(nMax);
        }

        virtual ~ParticleSet() {}
		/// Return the maximum number of particles allowed in memory.
        index_type nMax() const {return this->_nMax;}
        /// Returns the current number of particles (both active and inactive) in memory.
        index_type n() const {return this->_particles.size();}
        /// Returns the number of active particles (those that participate in the solution of a PDE)
        index_type nActive() const {return this->_nActive;}
        /// Returns the number of passive particles in memory (those that no longer participate in the solving of a PDE)
        index_type nPassive() const {return this->n() - this->_nActive;}
		/// The sum of the weights carried by all particles.  (Note: Passive particles should have their weight set to zero)
        scalar_type totalWeight() const;

		/// Return the set of all field names registered on these particles
        std::vector<std::string> fieldNames() const;
        
        inline bool isVertex(const index_type ind) const {return _particles[ind]->isVertex();}
        inline void isVertex(const index_type ind, const bool isvert) {_particles[ind]->isVertex(isvert);}

		/// Compute the integral of a scalar field defined on a ParticleSet
        virtual scalar_type scalarIntegral(const std::string& field_name) const;

		/// Initialize a ParticleSet from a VTK data file.
        void initFromParticleSetFile(const std::string& fname);

		/// Insert a new particle to the particle set
		/**
			New particles are added to the end of the array
		*/
        void insert(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0, const bool isvert=false);
        /// Insert a new particle to the particle set
		/**
			New particles are added to the end of the array
		*/
        void insert(const Vec<ndim>& xx, const scalar_type wght=0.0, const bool isvert=false);
        /// Move an existing particle to a new location in both physical and material space
        void move(const index_type ind, const Vec<ndim>& xx, const Vec<ndim>& aa);
        
        /// Return the ParticleFactory pointer associated with this ParticleSet
        inline ParticleFactory<ndim>* getFactory() const {return _factory.get();}

		/// Collect data abouta ParticleSet, return it in a string.
        virtual std::string infoString(const bool printAll=false) const;

		/// Return a Particle::infoString for each particle in the set.
        std::vector<std::string> particlesInfoStrings() const;

		/// Get the raw pointer to a particle in the set
        inline Particle<ndim>* getPtr(const index_type ind) const {return _particles[ind].get();}

		/// Initialize a scalar field on a particle set.
        void initScalarFieldFromFn(const std::string& field_name, const AnalyticFunction* fn);
        
        /// Intialize a vector field on a particle set.
        void initVectorFieldFromFn(const std::string& field_name, const AnalyticFunction* fn);
        
        /// Initialize a scalar field on a particle set.
        void initScalarFieldFromVector(const std::string& field_name, const std::vector<scalar_type>& vals);
        
        /// Intialize a vector field on a particle set.
        void initVectorFieldFromVectoryArray(const std::string& field_name, const std::vector<Vec<ndim>>& vals);
        
        /// Return a scalar field value associated with a particular particle
        std::vector<scalar_type> getScalarFieldValues(const std::string& field_name) const;
        
        inline scalar_type scalarVal(const index_type ind, const std::string field_name) const {
        	return _particles[ind]->getScalar(field_name);
        }
        /// Return a vector field value associated with a particular particle
        inline std::vector<scalar_type> vectorVal(const index_type ind, const std::string field_name) const {
        	return _particles[ind]->getVector(field_name);
        }
        
        /// Return the names of all registered fields (both scalar and vector).
        inline std::vector<std::string> getFieldNames() const {return _particles[0]->fieldNames();}
        /// Return the names of all scalar fields registered.
        inline std::vector<std::string> getScalarFieldNames() const {return _particles[0]->scalarFieldNames();}
        /// Return the names of all registered vector fields.
        inline std::vector<std::string> getVectorFieldNames() const {return _particles[0]->vectorFieldNames();}
		/// Return the physical coordinate of the particle at index ind
        inline Vec<ndim> physCrd(const index_type ind) const {return _particles[ind]->physCrd();}
        /// Return the Lagrangian coordinate of the particle at index ind
        inline Vec<ndim> lagCrd(const index_type ind) const {return _particles[ind]->lagCrd();}
        /// Return the weight of the particle at index ind
        inline scalar_type weight(const index_type ind) const {return _particles[ind]->weight();}
        /// Set the weight of the particle at index ind
        inline void setWeight(const index_type ind, const scalar_type val) { _particles[ind]->setWeight(val);}

#ifdef HAVE_VTK
		/// Convert particles to a VtkPoints instance
		vtkSmartPointer<vtkPoints> toVtkPoints(const bool useFieldForHeight = false, const std::string scalarFieldName="") const;
		/// Convert particle field data to a VtkPointData instance.
		vtkSmartPointer<vtkPointData> fieldsToVtkPointData() const;
		/// Convert particle field data to a VtkCellData instance.
		vtkSmartPointer<vtkCellData> fieldsToVtkCellData() const;
#endif

#ifdef HAVE_KOKKOS
		typedef Kokkos::View<scalar_type*[ndim]> vec_view_type;
		typedef Kokkos::View<scalar_type*> scalar_view_type;
		
		vec_view_type packPhysCoords() const;
		void unpackPhysCoords(const vec_view_type& cview);

// 		struct InitCoordView {
// 			vec_view_type cview;
// 			InitCoordView(vec_view_type cv) : cview(cv) {}
// 			KOKKOS_INLINE_FUNCTION
// 			void operator() (const index_type i) const {
// 				const Vec<ndim> pc 
// 				cview(i, 0) = 
// 			}
// 		};
		
		Kokkos::View<scalar_type*> packScalarField(const std::string& name) const;
		void unpackScalarField(const Kokkos::View<scalar_type*>& fview);
		
		Kokkos::View<scalar_type*[ndim]> packVectorField(const std::string& name) const;
		void unpackVectorField(const Kokkos::View<scalar_type*[ndim]>& fview);
#endif

    protected:
        index_type _nMax;
        index_type _nActive;
        std::shared_ptr<ParticleFactory<ndim>> _factory;
        std::vector<std::unique_ptr<Particle<ndim>>> _particles;
        
};
}
}
#endif
