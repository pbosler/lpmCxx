#include <iostream>
#include <sstream>
#include <fstream>
#include <exception>
#include "LpmAosParticleSet.hpp"

namespace Lpm {

template <int ndim> std::string ParticleSet<ndim>::infoString() const {
    std::ostringstream oss;
    oss << "particle set info:" << std::endl;
    oss << "\tnMax = " << _nMax << std::endl;
    oss << "\tn = " << _particles.size() << std::endl;
    oss << "\tnActive = " << _nActive << std::endl;
    oss << "\tfields:" << std::endl;
    std::vector<std::string> fields = fieldNames();
    for (int i=0; i<fields.size(); ++i)
        oss << "\t\t" << fields[i] << std::endl;
    oss << "\ttotalArea = " << totalArea() << std::endl;
    oss << "\ttotalVolume = " << totalVolume() << std::endl;
    return oss.str();
}

template <int ndim> std::vector<std::string> ParticleSet<ndim>::particlesInfoStrings() const {
    std::vector<std::string> result;
    for (index_type i=0; i<_particles.size(); ++i)
        result.push_back(_particles[i]->infoString());
    return result;
}

template <int ndim> scalar_type ParticleSet<ndim>::totalArea() const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_particles.size(); ++i)
        result += _particles[i]->area();
    return result;
}

template <int ndim> scalar_type ParticleSet<ndim>::totalVolume() const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_particles.size(); ++i)
        result += _particles[i]->volume();
    return result;
}

template <int ndim> std::vector<std::string> ParticleSet<ndim>::fieldNames() const {
    return _particles[0]->fieldNames();
}

template <int ndim> scalar_type ParticleSet<ndim>::scalarIntegral(const std::string& field_name) const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_particles.size(); ++i)
        result += _particles[i]->getScalar(field_name) * _particles[i]->area();
    return result;
}

template <int ndim> void ParticleSet<ndim>::initScalarFieldFromFn(const std::string& field_name,
    const AnalyticFunction* fn) {
    for (index_type i=0; i<_particles.size(); ++i) {
        _particles[i]->setScalar(field_name, fn->evaluateScalar(_particles[i]->physCrd()));
    }
}

template <int ndim> void ParticleSet<ndim>::initVectorFieldFromFn(const std::string& field_name,
    const AnalyticFunction* fn) {
    for (index_type i=0; i<_particles.size(); ++i) {
        _particles[i]->setVector(field_name, fn->evaluateVector(_particles[i]->physCrd()));
    }
}

template <int ndim> void ParticleSet<ndim>::insert(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type area,
    const scalar_type vol) {
    if (_particles.size() + 1 == _nMax) {
        throw std::out_of_range("ParticleSet nmax exceeded.");
    }
    _particles.push_back(_factory->createParticle(xx, aa, area, vol));
}

template <int ndim> void ParticleSet<ndim>::initFromParticleSetFile(const std::string& fname) {
    std::string fullfname(LPM_PARTICLE_SET_DIR);
    fullfname += "/" + fname;
    std::ifstream file(fullfname);
    if (!file.is_open()) {
        throw std::ios_base::failure("initFromParticleSetFile ERROR: file " + fullfname + " not found.");
    }

    index_type lineNumber = 0;
    std::string line;
    index_type nParticles = -1;
    bool areaFound = false;
    bool volumeFound = false;

    scalar_type xx, yy, zz, weight;
    index_type particleID = 0;
    while (std::getline(file, line)) {
        lineNumber += 1;
        if (lineNumber == 5) {
            nParticles = std::stol(line.substr(6));
            std::cout << "\t...found " << nParticles << " particles in file " << fullfname << std::endl;
            _nMax = nParticles;
            _particles.reserve(_nMax);
            _nActive = nParticles;
        }
        else if (nParticles > 0 && lineNumber > 5 && lineNumber <= 5 +nParticles ) {
            std::istringstream iss(line);
            _particles.push_back(_factory->createParticle());
            if (ndim == 2) {
                if (! (iss >> xx >> yy)) {
                    std::ostringstream oss;
                    oss << "initFromParticleSetFile read ERROR at line " << lineNumber;
                    oss << std::endl << "\t" << "nParticles = " << nParticles << std::endl;
                    oss << "\t" << line << std::endl;
                    throw std::ios_base::failure(oss.str());
                }
                _particles[particleID++]->init(Vec<ndim>(xx, yy));

            }
            else if (ndim == 3) {
                if (! (iss >> xx >> yy >> zz)) {
                    std::ostringstream oss;
                    oss << "initFromParticleSetFile read ERROR at line " << lineNumber;
                    oss << std::endl << "\t" << "nParticles = " << nParticles << std::endl;
                    oss << "\t" << line << std::endl;
                    throw std::ios_base::failure(oss.str());
                }
                _particles[particleID++]->init(Vec<ndim>(xx,yy,zz));
            }
        }
        else if (nParticles > 0 && lineNumber == 6 + nParticles) {
            particleID = 0;
        }
        else if (nParticles > 0 && lineNumber == 7 + nParticles) {
            if (line.find("area") != std::string::npos) {
                areaFound = true;
            }
            else if (line.find("volume") != std::string::npos) {
                volumeFound = true;
            }
        }
        else if (nParticles > 0 && lineNumber > 8 + nParticles) {
            std::istringstream iss(line);
            if (! (iss >> weight)) {
                std::ostringstream oss;
                oss << "initFromParticleSetFile read ERROR at line " << lineNumber;
                oss << std::endl << "\t" << "nParticles = " << nParticles << std::endl;
                oss << "\t" << line << std::endl;
                throw std::ios_base::failure(oss.str());
            }
            if (areaFound) {
                _particles[particleID++]->setArea(weight);
            }
            else if (volumeFound) {
                _particles[particleID++]->setVolume(weight);
            }
        }
    }

    file.close();
}

template class ParticleSet<3>;
template class ParticleSet<2>;
}
