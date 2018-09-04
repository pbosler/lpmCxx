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

namespace Lpm {
namespace Aos {

template <int ndim> class PolyMesh2d;

template <int ndim> class ParticleSet {
    public:
        ParticleSet() : _factory(nullptr), _nMax(0), _nActive(0) {}
    
        ParticleSet(const std::shared_ptr<ParticleFactory<ndim>> factory, const index_type nMax) : 
            _factory(factory), _nMax(nMax), _nActive(0) {
            _particles.reserve(nMax);
        }

        virtual ~ParticleSet() {}

        index_type nMax() const {return this->_nMax;}
        index_type n() const {return this->_particles.size();}
        index_type nActive() const {return this->_nActive;}
        index_type nPassive() const {return this->n() - this->_nActive;}

        scalar_type totalLength() const;
        scalar_type totalVolume() const;
        scalar_type totalArea() const;
        
//         inline bool particlesHaveLength() const {return _particles[0]->_hasLength;}
//         inline bool particlesHaveArea() const {return _particles[0]->_hasArea;}
//         inline bool particlesHaveVolume() const {return _particles[0]->_haveVolume;}

        std::vector<std::string> fieldNames() const;

        virtual scalar_type scalarIntegral(const std::string& field_name) const;

        void initFromParticleSetFile(const std::string& fname);

        void insert(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type len = 0., const scalar_type area=0.0, const scalar_type vol=0.0);
        void insert(const Vec<ndim>& xx, const scalar_type len = 0., const scalar_type area = 0.0, const scalar_type vol = 0.0);
        void move(const index_type ind, const Vec<ndim>& xx, const Vec<ndim>& aa);

        virtual std::string infoString() const;

        std::vector<std::string> particlesInfoStrings() const;

        inline Particle<ndim>* getPtr(const index_type ind) const {return _particles[ind].get();}

        void initScalarFieldFromFn(const std::string& field_name, const AnalyticFunction* fn);
        void initVectorFieldFromFn(const std::string& field_name, const AnalyticFunction* fn);
        
        void initScalarFieldFromVector(const std::string& field_name, const std::vector<scalar_type>& vals);
        void initVectorFieldFromVectoryArray(const std::string& field_name, const std::vector<Vec<ndim>>& vals);
        
        std::vector<scalar_type> getScalarFieldValues(const std::string& field_name) const;
        inline scalar_type scalarVal(const index_type ind, const std::string field_name) const {
        	return _particles[ind]->getScalar(field_name);
        }
        inline std::vector<scalar_type> vectorVal(const index_type ind, const std::string field_name) const {
        	return _particles[ind]->getVector(field_name);
        }
        
        inline std::vector<std::string> getFieldNames() const {return _particles[0]->fieldNames();}
        inline std::vector<std::string> getScalarFieldNames() const {return _particles[0]->scalarFieldNames();}
        inline std::vector<std::string> getVectorFieldNames() const {return _particles[0]->vectorFieldNames();}

        inline Vec<ndim> physCrd(const index_type ind) const {return _particles[ind]->physCrd();}
        inline Vec<ndim> lagCrd(const index_type ind) const {return _particles[ind]->lagCrd();}

#ifdef HAVE_VTK
		vtkSmartPointer<vtkPoints> toVtkPoints(const bool useFieldForHeight = false, const std::string scalarFieldName="") const;
		vtkSmartPointer<vtkPointData> fieldsToVtkPointData() const;
		vtkSmartPointer<vtkCellData> fieldsToVtkCellData() const;
#endif

		friend class PolyMesh2d<ndim>;
    protected:
        index_type _nMax;
        index_type _nActive;
        std::shared_ptr<ParticleFactory<ndim>> _factory;
        std::vector<std::unique_ptr<Particle<ndim>>> _particles;
        
};
}
}
#endif
