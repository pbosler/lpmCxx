#include <iostream>
#include <sstream>
#include <fstream>
#include <exception>
#include "LpmAosParticleSet.hpp"
#ifdef HAVE_VTK
#include "vtkDoubleArray.h"
#endif

namespace Lpm {
namespace Aos {

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
    oss << "\ttotal "<< _particles[0]->weightName() << " = " << totalWeight() << std::endl;
    return oss.str();
}

template <int ndim> std::vector<std::string> ParticleSet<ndim>::particlesInfoStrings() const {
    std::vector<std::string> result;
    for (index_type i=0; i<_particles.size(); ++i)
        result.push_back(_particles[i]->infoString());
    return result;
}

template <int ndim> scalar_type ParticleSet<ndim>::totalWeight() const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_particles.size(); ++i) {
        result += _particles[i]->weight();
    }
    return result;
}

template <int ndim> std::vector<std::string> ParticleSet<ndim>::fieldNames() const {
    return (!_particles.empty() ? _particles[0]->fieldNames() : std::vector<std::string>(1, "null"));
}

template <int ndim> scalar_type ParticleSet<ndim>::scalarIntegral(const std::string& field_name) const {
    scalar_type result = 0.0;
    for (index_type i=0; i<_particles.size(); ++i)
        result += _particles[i]->getScalar(field_name) * _particles[i]->weight();
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

template <int ndim> void ParticleSet<ndim>::initScalarFieldFromVector(const std::string& field_name, const std::vector<scalar_type>& vals) {
    if (_nActive > vals.size()) {
        throw std::out_of_range("ParticleSet::initScalarFieldFromVector ERROR : size mismatch.");
    }
    for (index_type i=0; i<_nActive; ++i) {
        _particles[i]->setScalar(field_name, vals[i]);
    }
}

template <int ndim> void ParticleSet<ndim>::initVectorFieldFromVectoryArray(const std::string& field_name, const std::vector<Vec<ndim>>& vals) {
    if (_nActive > vals.size()) {
        throw std::out_of_range("ParticleSet::initVectorFieldFromVectorArray ERROR : size mismatch.");
    }
    for (index_type i=0; i<_nActive; ++i) {
        _particles[i]->setVector(field_name, vals[i]);
    }
}

template <int ndim> void ParticleSet<ndim>::insert(const Vec<ndim>& xx, const Vec<ndim>& aa, const scalar_type wgt) {
    if (_particles.size() + 1 > _nMax) {
        throw std::out_of_range("ParticleSet nmax exceeded.");
    }
    _particles.push_back(_factory->createParticle(xx, aa, wgt));
    _nActive += 1;
}

template <int ndim> void ParticleSet<ndim>::insert(const Vec<ndim>& xx, const scalar_type wgt) {
    if (_particles.size() + 1 > _nMax) {
        throw std::out_of_range("ParticleSet nmax exceeded.");
    }
    _particles.push_back(_factory->createParticle(xx, wgt));    
    _nActive += 1;
}

template <int ndim> void ParticleSet<ndim>::move(const index_type ind, const Vec<ndim>& xx, const Vec<ndim>& aa) {
    _particles[ind]->move(xx, aa);
}

template <int ndim> std::vector<scalar_type> ParticleSet<ndim>::getScalarFieldValues(const std::string& field_name) const {
    std::vector<scalar_type> result(_nActive);
    for (index_type i=0; i<_nActive; ++i) {
        result[i] = _particles[i]->getScalar(field_name);
    }
    return result;
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
            //_nMax = nParticles;
            _particles.reserve(nParticles);
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
                _particles[particleID++]->setWeight(weight);
            }
            else if (volumeFound) {
                _particles[particleID++]->setWeight(weight);
            }
        }
    }

    file.close();
}

#ifdef HAVE_VTK
template <int ndim> vtkSmartPointer<vtkPoints> ParticleSet<ndim>::toVtkPoints(const bool useFieldForHeight, 
	const std::string scalarFieldName) const {
		vtkSmartPointer<vtkPoints> result = vtkSmartPointer<vtkPoints>::New();
		switch (ndim) {
			case (2) : {
				if (useFieldForHeight) {
					for (index_type i=0; i<_nActive; ++i) {
						const Particle<ndim>* pptr = getPtr(i);
						const Vec<ndim> pos=pptr->physCrd();
						result->InsertPoint(i, pos.x[0], pos.x[1], pptr->getScalar(scalarFieldName));
					}
				}
				else {
					for (index_type i=0; i<_nActive; ++i) {
						const Vec<ndim> pos=getPtr(i)->physCrd();
						result->InsertPoint(i, pos.x[0], pos.x[1], 0.0);
					}
				}
				break;
			}
			case (3) : {
				for (index_type i=0; i<_nActive; ++i) {
					const Vec<ndim> pos=getPtr(i)->physCrd();
					result->InsertPoint(i, pos.x[0], pos.x[1], pos.x[2]);
				}
			}
		}
		return result;
	}
	
template <int ndim>	vtkSmartPointer<vtkPointData> ParticleSet<ndim>::fieldsToVtkPointData() const{
	vtkSmartPointer<vtkPointData> result = vtkSmartPointer<vtkPointData>::New();
	// add geometric data
	vtkSmartPointer<vtkDoubleArray> wgt = vtkSmartPointer<vtkDoubleArray>::New();
	wgt->SetName(_particles[0]->weightName().c_str());
	wgt->SetNumberOfComponents(1);
	wgt->SetNumberOfTuples(_nActive);
	for (index_type j=0; j<_nActive; ++j) {
		wgt->InsertTuple1(j, _particles[j]->weight());
	}
	result->AddArray(wgt);

	// collect field names
	const std::vector<std::string> sfields = getScalarFieldNames();
	const std::vector<std::string> vfields = getVectorFieldNames();
	// add field data
	for (int i=0; i<sfields.size(); ++i) {
	    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
	    data->SetName(sfields[i].c_str());
	    data->SetNumberOfComponents(1);
	    data->SetNumberOfTuples(_nActive);
	    for (index_type j=0; j<_nActive; ++j) {
	        data->InsertTuple1(j, _particles[j]->getScalar(sfields[i]));
	    }
	    result->AddArray(data);
	}   
	for (int i=0; i<vfields.size(); ++i) {
	    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
	    data->SetName(vfields[i].c_str());
	    data->SetNumberOfComponents(ndim);
	    data->SetNumberOfTuples(_nActive);
	    switch (ndim) {
	        case (2) : {
	            for (index_type j=0; j<_nActive; ++j) {
	                const std::vector<scalar_type> val = _particles[j]->getVector(vfields[i]);
	                data->InsertTuple2(j, val[0], val[1]);
	            }
	            break;
	        }
	        case (3) : {
	            for (index_type j=0; j<_nActive; ++j) {
	                const std::vector<scalar_type> val = _particles[j]->getVector(vfields[i]);
	                data->InsertTuple3(j, val[0], val[1], val[2]);
	            }
	            break;
	        }
	    }
	    result->AddArray(data);
	}
	return result;
}

template <int ndim> vtkSmartPointer<vtkCellData> ParticleSet<ndim>::fieldsToVtkCellData() const {
    vtkSmartPointer<vtkCellData> result = vtkSmartPointer<vtkCellData>::New();
    // Add geometric quantities
	vtkSmartPointer<vtkDoubleArray> wgt = vtkSmartPointer<vtkDoubleArray>::New();
	wgt->SetName(_particles[0]->weightName().c_str());
	wgt->SetNumberOfComponents(1);
	wgt->SetNumberOfTuples(_nActive);
	for (index_type j=0; j<_nActive; ++j) {
		wgt->InsertTuple1(j, _particles[j]->weight());
	}
	result->AddArray(wgt);
    // collect field names
    const std::vector<std::string> sfields = getScalarFieldNames();
	const std::vector<std::string> vfields = getVectorFieldNames();
	// add field data
	for (int i=0; i<sfields.size(); ++i) {
	    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
	    data->SetName(sfields[i].c_str());
	    data->SetNumberOfComponents(1);
	    data->SetNumberOfTuples(_nActive);
	    for (index_type j=0; j<_nActive; ++j) {
	        data->InsertTuple1(j, _particles[j]->getScalar(sfields[i]));
	    }
	    result->AddArray(data);
	}   
	for (int i=0; i<vfields.size(); ++i) {
	    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
	    data->SetName(vfields[i].c_str());
	    data->SetNumberOfComponents(ndim);
	    data->SetNumberOfTuples(_nActive);
	    switch (ndim) {
	        case (2) : {
	            for (index_type j=0; j<_nActive; ++j) {
	                const std::vector<scalar_type> val = _particles[j]->getVector(vfields[i]);
	                data->InsertTuple2(j, val[0], val[1]);
	            }
	            break;
	        }
	        case (3) : {
	            for (index_type j=0; j<_nActive; ++j) {
	                const std::vector<scalar_type> val = _particles[j]->getVector(vfields[i]);
	                data->InsertTuple3(j, val[0], val[1], val[2]);
	            }
	            break;
	        }
	    }
	    result->AddArray(data);
	}
    return result;
}

#endif

template class ParticleSet<1>;
template class ParticleSet<2>;
template class ParticleSet<3>;
}
}