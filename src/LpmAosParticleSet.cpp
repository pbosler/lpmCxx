#include <iostream>
#include <sstream>
#include <fstream>
#include <exception>
#include "LpmAosParticleSet.hpp"

namespace Lpm {

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
                    nParticles = std::stol(line.substr(6,2));
                }
                else if (lineNumber > 5 && lineNumber < 6 +nParticles ) {
                    std::istringstream iss(line);
                    _particles.push_back(_factory->createParticle());
                    if (ndim == 2) {
                        if (! (iss >> xx >> yy)) {
                            std::ostringstream oss;
                            oss << "initFromParticleSetFile read ERROR at line " << lineNumber;
                            throw std::ios_base::failure(oss.str());
                        }
                        _particles[particleID++]->init(Vec<ndim>(xx, yy));
                
                    }
                    else if (ndim == 3) {
                        if (! (iss >> xx >> yy >> zz)) {
                            std::ostringstream oss;
                            oss << "initFromParticleSetFile read ERROR at line " << lineNumber;
                            throw std::ios_base::failure(oss.str());
                        }
                        _particles[particleID++]->init(Vec<ndim>(xx,yy,zz));
                    }
                }
                else if (lineNumber == 6 + nParticles) {
                    particleID = 0;
                }
                else if (lineNumber == 7 + nParticles) {
                    if (line.find("area") != std::string::npos) {
                        areaFound = true;
                    }
                    else if (line.find("volume") != std::string::npos) {
                        volumeFound = true;
                    }
                }
                else if (lineNumber > 8 + nParticles) {
                    std::istringstream iss(line);
                    if (! (iss >> weight)) {
                        std::ostringstream oss;
                        oss << "initFromParticleSetFile read ERROR at line " << lineNumber;
                        throw std::ios_base::failure(oss.str());
                    }
                    if (areaFound) {
                        _particles[particleID]->setArea(weight);
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
