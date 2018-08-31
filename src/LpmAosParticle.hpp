#ifndef _LPM_AOS_PARTICLE_HPP
#define _LPM_AOS_PARTICLE_HPP

#include "LpmAosTypes.hpp"
#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include <string>
#include <sstream>
#include <map>
#include <array>

namespace Lpm {

template <int ndim=3> class Particle {
    typedef std::array<scalar_type, ndim>  vfield_type;

    protected:
        Vec<ndim> _physCrd;
        Vec<ndim> _lagCrd;
        vfield_type _velocity;
        scalar_type _length;
        scalar_type _area;
        scalar_type _volume;

        std::map<std::string, scalar_type> _sFields;
        std::map<std::string, vfield_type> _vFields;

    public:

        Particle() : _physCrd(), _lagCrd(), _area(0.0), _volume(0.0), _velocity() {}
        Particle(const Vec<ndim>& pos, const scalar_type len = 0.0, const scalar_type aa=0.0, const scalar_type vv=0.0) : _physCrd(pos), _lagCrd(pos),
             _length(len), _area(aa), _volume(vv), _velocity() {}
        Particle(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type len=0.0, const scalar_type ar=0.0, const scalar_type vol=0.0) :
            _physCrd(xx), _lagCrd(aa), _area(ar), _volume(vol), _velocity() {}

        void setLength(const scalar_type l) {_length=l;}
        void setArea(const scalar_type a) {_area = a;}
        void setVolume(const scalar_type v) {_volume = v;}

        inline Vec<ndim> physCrd() const {return _physCrd;}
        inline Vec<ndim> lagCrd() const {return _lagCrd;}

        inline scalar_type length() const {return _length;}
        inline scalar_type area() const {return _area;}
        inline scalar_type volume() const {return _volume;}

        void registerScalarField(const std::string& field_name) {
            this->_sFields.emplace(field_name, 0.0);
        }

        void registerVectorField(const std::string& field_name) {
            this->_vFields.emplace(field_name, vfield_type());
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
            _vFields.at(field_name) = vec.toArray();
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

        virtual void init(const Vec<ndim>& initCrd, const scalar_type len=0.0, const scalar_type aa=0.0, const scalar_type vv=0.0) {
            _physCrd = initCrd;
            _lagCrd = initCrd;
            _length = len;
            _area = aa;
            _volume = vv;
        };
        
        virtual void move(const Vec<ndim>& xx, const Vec<ndim>& aa) {
            _physCrd = xx;
            _lagCrd = aa;
        }

        virtual void init(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type len=0.0, const scalar_type ar=0.0, const scalar_type vv=0.0) {
            _physCrd = xx;
            _length = len;
            _lagCrd = aa;
            _area = ar;
            _volume = vv;
        };

        virtual ~Particle() {}

        void setVel(const scalar_type* vv) {
            for (int i=0; i<ndim; ++i)
                this->velocity[i] = vv[i];
        }

        void setVel(const std::vector<scalar_type>& vv) {
            for (int i=0; i<ndim; ++i)
                this->_velocity[i] = vv[i];
        }

        void setVel(const scalar_type u, const scalar_type v) {
            this->_velocity[0] = u;
            this->_velocity[1] = v;
        }

        void setVel(const scalar_type u, const scalar_type v, const scalar_type w) {
            this->_velocity[0] = u;
            this->_velocity[1] = v;
            this->_velocity[2] = w;
        }

        std::string infoString() const {
            std::ostringstream ss;
            ss << "particle info:" << std::endl;
            ss << "\tphysCrd = " << this->_physCrd << std::endl;
            ss << "\tlagCrd = " << this->_lagCrd << std::endl;
            ss << "\tvelocity = [ ";
            for (int i=0; i<ndim; ++i)
                ss << this->_velocity[i] << " ";
            ss << "]" << std::endl;
            ss << "\tlength = " << this->_length << std::endl;
            ss << "\tarea = " << this->_area << std::endl;
            ss << "\tvolume = " << this->_volume << std::endl;
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

#endif
