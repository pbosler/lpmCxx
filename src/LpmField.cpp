#include "LpmField.h"
#include <sstream>
#include <algorithm>
#include <exception>
#include <numeric>
#include <cmath>

namespace Lpm {

std::unique_ptr<Logger> Field::log(new Logger(OutputMessage::debugPriority, "Field_log"));

Field::Field(const index_type nMax, const int nDim, const std::string name, const std::string units) : 
    _nMax(nMax), _nDim(nDim), _name(name), _units(units) {
    switch (nDim) {
        case (1) : {
            comp0.reserve(nMax);
            break;
        }
        case (2) : {
            comp0.reserve(nMax);
            comp1.reserve(nMax);
            break;
        }
        case (3) : {
            comp0.reserve(nMax);
            comp1.reserve(nMax);
            comp2.reserve(nMax);
            break;
        }
        default : {
            OutputMessage errMsg("ERROR invalid dimensions (nDim)", OutputMessage::errorPriority, "Field::Field");
            log->logMessage(errMsg);
            throw std::runtime_error("invalid input");
        }
    }
}

scalar_type Field::maxScalarVal() const {
    return *std::max_element(comp0.begin(), comp0.end());
}

scalar_type Field::minScalarVal() const {
    return *std::min_element(comp0.begin(), comp0.end());
}

std::string Field::listFieldValues() const {
    std::stringstream ss;
    switch (_nDim) {
        case (1) : {
            for (index_type i = 0; i < comp0.size(); ++i) {
                ss << comp0[i] << std::endl;
            }
            break;
        }
        case (2) : {
            for (index_type i = 0; i < comp0.size(); ++i) {
                const XyzVector vec = XyzVector(comp0[i], comp1[i]);
                ss << vec << std::endl;
            }
            break;
        }
        case (3) : {
            for (index_type i = 0; i < comp0.size(); ++i) {
                const XyzVector vec = XyzVector(comp0[i], comp1[i], comp2[i]);
                ss << vec << std::endl;
            }
            break;
        }
    }
    return ss.str();
}

void Field::abs() {
    switch (_nDim) {
        case (1) : {
            for (index_type i = 0; i < n(); ++i) {
                comp0[i] = std::abs(comp0[i]);
            }
            break;
        }
        case (2) : { 
            for (index_type i = 0; i < n(); ++i) {
                comp0[i] = std::abs(comp0[i]);
                comp1[i] = std::abs(comp1[i]);
            }
            break;
        }
        case (3) : {
            for (index_type i = 0; i < n(); ++i) {
                comp0[i] = std::abs(comp0[i]);
                comp1[i] = std::abs(comp1[i]);
                comp2[i] = std::abs(comp2[i]);
            }
            break;
        }
    }
}


void Field::scale(const scalar_type multiplier) {
    switch (_nDim) {
        case (1) : {
            for (index_type i = 0; i < n(); ++i) {
                comp0[i] *= multiplier;
            }
            break;
        }
        case (2) : { 
            for (index_type i = 0; i < n(); ++i) {
                comp0[i] *= multiplier;
                comp1[i] *= multiplier;
            }
            break;
        }
        case (3) : {
            for (index_type i = 0; i < n(); ++i) {
                comp0[i] *= multiplier;
                comp1[i] *= multiplier;
                comp2[i] *= multiplier;
            }
            break;
        }
    }
}

void Field::update(const scalar_type a, const std::shared_ptr<Field>& other, 
                   const scalar_type b, const std::shared_ptr<Field>& other1) {
    if (_nDim != other->nDim()) {
        std::stringstream ss;
        ss << _name << " field " << " does not have same dimension as " << other->_name << "; cannot update.";
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Field::update");
        log->logMessage(errMsg);
        throw std::runtime_error("field dimension mismatch");
    }
    if ( n() < other->n() ) {
        std::stringstream ss;
        ss << _name << " field " << " has insufficient memory to update with field " << other->_name;
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Field::update");
        log->logMessage(errMsg);
        throw std::runtime_error("field memory mismatch");
    }
    if ( other1 && b != 0.0) {
        if (_nDim != other1->nDim()) {
            std::stringstream ss;
            ss << _name << " field " << " does not have same dimension as " << other1->_name << "; cannot update.";
            OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Field::update");
            log->logMessage(errMsg);
            throw std::runtime_error("field dimension mismatch");
        }
        if ( n() < other1->n() ) {
            std::stringstream ss;
            ss << _name << " field " << " has insufficient memory to update with field " << other1->_name;
            OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Field::update");
            log->logMessage(errMsg);
            throw std::runtime_error("field memory mismatch");
        }
    }
    
    if (other1 && b != 0.0) {
        switch (_nDim) {
            case (1) : {
                for (index_type i = 0; i < n(); ++i) {
                    comp0[i] = a * other->comp0[i] + b * other1->comp0[i];
                }
                break;
            }
            case (2) : { 
                for (index_type i = 0; i < n(); ++i) {
                    comp0[i] = a * other->comp0[i] + b * other1->comp0[i];
                    comp1[i] = a * other->comp0[i] + b * other1->comp1[i];
                }
                break;
            }
            case (3) : {
                for (index_type i = 0; i < n(); ++i) {
                    comp0[i] = a * other->comp0[i] + b * other1->comp0[i];
                    comp1[i] = a * other->comp1[i] + b * other1->comp1[i];
                    comp2[i] = a * other->comp2[i] + b * other1->comp2[i];
                }
                break;
            }
        }
    }
    else {
        switch (_nDim) {
            case (1) : {
                for (index_type i = 0; i < n(); ++i) {
                    comp0[i] = a * other->comp0[i];
                }
                break;
            }
            case (2) : { 
                for (index_type i = 0; i < n(); ++i) {
                    comp0[i] = a * other->comp0[i];
                    comp1[i] = a * other->comp1[i];
                }
                break;
            }
            case (3) : {
                for (index_type i = 0; i < n(); ++i) {
                    comp0[i] = a * other->comp0[i];
                    comp1[i] = a * other->comp1[i];
                    comp2[i] = a * other->comp2[i];
                }
                break;
            }
        }
    }
}

void Field::insert(const scalar_type fx) {
    if (n() + 1 > _nMax) {
        std::stringstream ss;
        ss << "not enough memory to insert scalar " << fx;
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Field::insert(scalar)");
        log->logMessage(errMsg);
        throw std::bad_alloc();
    }
    comp0.push_back(fx);
}

void Field::insert(const scalar_type fx, const scalar_type fy, const scalar_type fz) {
    if (n() + 1 > _nMax) {
        std::stringstream ss;
        ss << "not enough memory to insert vector " << XyzVector(fx, fy, fz);
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Field::insert(3 scalars)");
        log->logMessage(errMsg);
        throw std::bad_alloc();
    }
    comp0.push_back(fx);
    comp1.push_back(fy);
    if (_nDim == 3)
        comp2.push_back(fz);
}

void Field::insert(const XyzVector& vec) {
    if (n() + 1 > _nMax) {
        std::stringstream ss;
        ss << "not enough memory to insert vector " << vec;
        OutputMessage errMsg(ss.str(), OutputMessage::errorPriority, "Field::insert(XyzVector)");
        log->logMessage(errMsg);
        throw std::bad_alloc();
    }
    comp0.push_back(vec.x);
    comp1.push_back(vec.y);
    if (_nDim == 3)
        comp2.push_back(vec.z);
}

void Field::replace(const index_type ind, const scalar_type fx) {
    comp0[ind] = fx;
}

void Field::replace(const index_type ind, const scalar_type fx, const scalar_type fy, const scalar_type fz) {
    comp0[ind] = fx;
    comp1[ind] = fy;
    if (_nDim == 3)
        comp2[ind] = fz;
}

void Field::replace(const index_type ind, const XyzVector& vec) {
    comp0[ind] = vec.x;
    comp1[ind] = vec.y;
    if (_nDim == 3)
        comp2[ind] = vec.z;
}

void Field::initializeToConstant(const Coords* crds, const scalar_type val) {
    switch (_nDim) {
        case (1) : {
            comp0 = std::vector<scalar_type>(crds->n(), val);
            break;
        }
        case (2) : {
            comp0 = std::vector<scalar_type>(crds->n(), val);
            comp1 = std::vector<scalar_type>(crds->n(), val);
            break;
        }
        case (3) : {
            comp0 = std::vector<scalar_type>(crds->n(), val);
            comp1 = std::vector<scalar_type>(crds->n(), val);
            comp2 = std::vector<scalar_type>(crds->n(), val);
            break;
        }
    }
}

 
void Field::initializeToConstant(const index_type nn, const scalar_type val) {
    switch (_nDim) {
        case (1) : {
            comp0 = std::vector<scalar_type>(nn, val);
            break;
        }
        case (2) : {
            comp0 = std::vector<scalar_type>(nn, val);
            comp1 = std::vector<scalar_type>(nn, val);
            break;
        }
        case (3) : {
            comp0 = std::vector<scalar_type>(nn, val);
            comp1 = std::vector<scalar_type>(nn, val);
            comp2 = std::vector<scalar_type>(nn, val);
            break;
        }
    }
} 

void Field::initializeToScalarFunction(const Coords* crds, const AnalyticFunction* fn) {
    comp0.clear();
    for (index_type i = 0; i < crds->n(); ++i) {
        comp0.push_back(fn->evaluateScalar(crds->getVec(i)));
    }
}

void Field::initializeToVectorFunction(const Coords* crds, const AnalyticFunction* fn) {
    comp0.clear();
    comp1.clear();
    if (_nDim == 3)
        comp2.clear();
    for (index_type i = 0; i < crds->n(); ++i) {
        const XyzVector vecVal = fn->evaluateVector(crds->getVec(i));
        comp0.push_back(vecVal.x);
        comp1.push_back(vecVal.y);
        if (_nDim == 3)
            comp2.push_back(vecVal.z);
    }
}

}
