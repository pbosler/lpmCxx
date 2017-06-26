#include "LpmField.h"

namespace Lpm {

Field::Field(const index_type nMax, const int nDim, const std::string name, const std::string units) : 
    _nMax(nMax), _nDim(nDim), _name(name), _units(units) {
    if (!log) {
        log = std::unique_ptr<Logger>(new Logger(OutputMessage::debugPriority));
    }
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
            break;
        }
    }
}

void Field::insert(const scalar_type fx) {
    comp0.push_back(fx);
}

void Field::insert(const scalar_type fx, const scalar_type fy, const scalar_type fz) {
    comp0.push_back(fx);
    comp1.push_back(fy);
    if (_nDim == 3)
        comp2.push_back(fz);
}

void Field::insert(const XyzVector& vec) {
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

void Field::initializeToConstant(const scalar_type val) {
    switch (_nDim) {
        case (1) : {
            for (index_type i = 0; i < comp0.size(); ++i) {
                comp0[i] = val;
            }
            break;
        }
        case (2) : {
            for (index_type i = 0; i < comp0.size(); ++i) {
                comp0[i] = val;
                comp1[i] = val;
            }
            break;
        }
        case (3) : {
            for (index_type i = 0; i < comp0.size(); ++i) {
                comp0[i] = val;
                comp1[i] = val;
                comp2[i] = val;
            }
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
