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
    Then its values may be initialized.  They may then be modified by other classes (solvers, e.g.)
*/
template <int ndim=3> class Particle {
    typedef std::vector<scalar_type>  vfield_type;

    protected:
        Vec<ndim> _physCrd;
        Vec<ndim> _lagCrd;
        scalar_type _weight;
        std::string _wgt_name;

        std::map<std::string, scalar_type> _sFields;
        std::map<std::string, vfield_type> _vFields;

    public:
        //friend class ParticleSet<ndim>;

        Particle(const std::string& wname="nullweight") : _physCrd(), _lagCrd(), _weight(0.0), _wgt_name(wname) {}
        Particle(const Vec<ndim>& pos, const scalar_type wgt = 0.0, const std::string wname="nullweight") : _physCrd(pos), _lagCrd(pos), _weight(wgt), _wgt_name(wname) {}
        Particle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0, const std::string wname="nullweight") : _physCrd(xx), _lagCrd(aa), _weight(wgt), _wgt_name(wname){}

        

        inline Vec<ndim> physCrd() const {return _physCrd;}
        inline Vec<ndim> lagCrd() const {return _lagCrd;}

        inline void setWeight(const scalar_type wgt) {_weight = wgt;}
        inline scalar_type weight() const {return _weight;}
        inline std::string weightName() const {return _wgt_name;}

        void registerScalarField(const std::string& field_name) {
            this->_sFields.emplace(field_name, 0.0);
        }

        void registerVectorField(const std::string& field_name) {
            this->_vFields.emplace(field_name, vfield_type(ndim, 0.0));
        }
        
        std::vector<std::string> scalarFieldNames() const {
            std::vector<std::string> result;
            for (auto& sf : _sFields) {
                result.push_back(sf.first);
            }
            return result;
        }
        std::vector<std::string> vectorFieldNames() const {
            std::vector<std::string> result;
            for (auto& vf : _vFields) {
                result.push_back(vf.first);
            }
            return result;
        }

        std::vector<std::string> fieldNames() const {
            std::vector<std::string> result;
            for (auto& vf : _vFields)
                result.push_back(vf.first);
            for (auto& sf : _sFields)
                result.push_back(sf.first);
            return result;
        }

        inline void setScalar(const std::string& field_name, const scalar_type val) {
            _sFields.at(field_name) = val;
        }

        inline void setVector(const std::string& field_name, const Vec<ndim>& vec) {
            _vFields.at(field_name) = vec.toStdVec();
        }

        inline void setVector(const std::string& field_name, const vfield_type& arr) {
            _vFields[field_name] = arr;
        }

        scalar_type getScalar(const std::string& field_name) const {
            return _sFields.at(field_name);
        }

        vfield_type getVector(const std::string& field_name) const {
            return _vFields.at(field_name);
        }

        virtual void init(const Vec<ndim>& initCrd, const scalar_type wgt=0.0) {
            _physCrd = initCrd;
            _lagCrd = initCrd;
            _weight = wgt;
        };
        
        virtual void move(const Vec<ndim>& xx, const Vec<ndim>& aa) {
            _physCrd = xx;
            _lagCrd = aa;
        }

        virtual void init(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt=0.0, const std::string& wname="null") {
            _physCrd = xx;
            _lagCrd = aa;
            _weight = wgt;
            _wgt_name = wname;
        }

        virtual ~Particle() {}

        std::string infoString() const {
            std::ostringstream ss;
            ss << "particle info:" << std::endl;
            ss << "\tphysCrd = " << this->_physCrd << std::endl;
            ss << "\tlagCrd = " << this->_lagCrd << std::endl;
//             ss << "\tvelocity = [ ";
//             for (int i=0; i<ndim; ++i)
//                 ss << this->_velocity[i] << " ";
//             ss << "]" << std::endl;
            ss << "\tweight(" << _wgt_name << ") = " << this->_weight << std::endl;
            for (auto& sf : this->_sFields) {
                ss << "\t" << sf.first << " = " << sf.second << std::endl;
            }
            for (auto& vf : this->_vFields) {
                ss << "\t" << vf.first << " = [ ";
                for (int i=0; i<ndim; ++i)
                    ss << vf.second[i] << " ";
                ss << "]" << std::endl;
            }
            return ss.str();
        }
};

}
}
#endif
