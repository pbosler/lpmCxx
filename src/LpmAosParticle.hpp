#ifndef _LPM_AOS_PARTICLE_HPP
#define _LPM_AOS_PARTICLE_HPP

#include "LpmAosTypes.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include <string>
#include <sstream>
#include <map>
#include <vector>

namespace Lpm {
namespace Aos {

template <int ndim> class ParticleSet;

/// Basic discretization unit; A particle may represent any data distributed in R^d, where d \in {1,2,3}.
/*
    Data are carried by particles.  First, each data field must be "registered," as either a scalar or a vector.
    Then its values may be initialized.  They may then be modified by other classes (solvers, e.g.).
    
    Each particle has a coordinate in physical space and a coordinate in Lagrangian (or material) space.  
    Additionally, each particle carries a  "weight," e.g., area for a 2D problem, volume for a 3d problem.
    These weights may be generalized to mass densities and/or quadrature weights, etc.
*/
template <int ndim=3> class Particle {
    //typedef std::array<scalar_type, ndim>  vfield_type;
    typedef std::vector<scalar_type> vfield_type;

    protected:
        Vec<ndim> _physCrd;
        Vec<ndim> _lagCrd;
        scalar_type _weight;
        bool _is_vertex;

        std::map<std::string, scalar_type> _sFields;
        std::map<std::string, vfield_type> _vFields;

    public:
        virtual ~Particle() {}

        /// Constructor.  Initializes to zero/null.
        Particle(const bool isvert=false) : _physCrd(), _lagCrd(), _weight(0.0), _is_vertex(isvert) {}
        /// Constructor.  Physical position and weight only.
        Particle(const Vec<ndim>& pos, const scalar_type wgt = 0.0, const bool isvert=false) : 
            _physCrd(pos), _lagCrd(pos), _weight(wgt), _is_vertex(isvert) {}
        /// Constructor.  Physical and Lagrangian positions and weight.
        Particle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0,
            const bool isvert=false) : _physCrd(xx), _lagCrd(aa), _weight(wgt), _is_vertex(isvert) {}


        inline bool isVertex() const {return _is_vertex;}
        inline void isVertex(const bool isvert) {_is_vertex = isvert;}
        
        /// Return a particle's physcial coordinate
        inline Vec<ndim> physCrd() const {return _physCrd;}
        /// Return a particle's Lagrangian coordinate
        inline Vec<ndim> lagCrd() const {return _lagCrd;}
        
        /// Set a particle's geometric weight
        inline void setWeight(const scalar_type wgt) {_weight = wgt;}
        /// Return a particle's weight
        inline scalar_type weight() const {return _weight;}


        /// Register a new scalar field on this particle
        void registerScalarField(const std::string& field_name) {
            this->_sFields.emplace(field_name, 0.0);
        }

        /// Register a new vector field on this particle
        void registerVectorField(const std::string& field_name) {
            this->_vFields.emplace(field_name, vfield_type(ndim, 0.0));
        }
        
        /// Return the names of all scalar fields registered.
        std::vector<std::string> scalarFieldNames() const {
            std::vector<std::string> result;
            for (auto& sf : _sFields) {
                result.push_back(sf.first);
            }
            return result;
        }
        /// Return the names of all registered vector fields.
        std::vector<std::string> vectorFieldNames() const {
            std::vector<std::string> result;
            for (auto& vf : _vFields) {
                result.push_back(vf.first);
            }
            return result;
        }

        /// Return the names of all registered fields (both scalar and vector).
        std::vector<std::string> fieldNames() const {
            std::vector<std::string> result;
            for (auto& vf : _vFields)
                result.push_back(vf.first);
            for (auto& sf : _sFields)
                result.push_back(sf.first);
            return result;
        }

        /// Set the value of a scalar field at this particle
        inline void setScalar(const std::string& field_name, const scalar_type val) {
            _sFields.at(field_name) = val;
        }

        /// Set the value of a vector field at this particle
        inline void setVector(const std::string& field_name, const Vec<ndim>& vec) {
            _vFields.at(field_name) = vec.toStdVec();
        }

        /// Set the value of a vector field at this particle
        inline void setVector(const std::string& field_name, const vfield_type& arr) {
            _vFields[field_name] = arr;
        }
    
        /// Get the value of a scalar field at this particle
        scalar_type getScalar(const std::string& field_name) const {
            return _sFields.at(field_name);
        }
        
        /// Get the value of a vector field at this particle
        vfield_type getVector(const std::string& field_name) const {
            return _vFields.at(field_name);
        }

        /// Initialize a particle
        /**
            @deprecated
        */
        virtual void init(const Vec<ndim>& initCrd, const scalar_type wgt=0.0) {
            _physCrd = initCrd;
            _lagCrd = initCrd;
            _weight = wgt;
        };
        
        /// Move a particle to new coordinates
        virtual void move(const Vec<ndim>& xx, const Vec<ndim>& aa) {
            _physCrd = xx;
            _lagCrd = aa;
        }

        /// Initialize a particle
        /**
            @deprecated
        */
        virtual void init(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0) {
            _physCrd = xx;
            _lagCrd = aa;
            _weight = wgt;
        }

        

        /// Return a string with data about the particle.  
        std::string infoString() const {
            std::ostringstream ss;
            ss << "particle info:" << std::endl;
            ss << "\tphysCrd = " << this->_physCrd << std::endl;
            ss << "\tlagCrd = " << this->_lagCrd << std::endl;
            ss << "\tweight = " << this->_weight << std::endl;
            for (auto& sf : this->_sFields) {
                ss << "\t" << sf.first << " = " << sf.second << std::endl;
            }
            for (auto& vf : this->_vFields) {
                ss << "\t" << vf.first << " = [ ";
                for (int i=0; i<ndim; ++i)
                    ss << vf.second[i] << " ";
                ss << "]" << std::endl;
            }
            ss << "\tis vertex? " << (_is_vertex ? "yes" : "no") << std::endl;
            return ss.str();
        }
};

}
}
#endif
