#ifndef _LPM_FIELD_H_
#define _LPM_FIELD_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmLogger.h"
#include "LpmAnalyticFunctions.h"
#include "LpmCoords.h"
#include <vector>
#include <string>
#include <memory>

namespace Lpm {

class Field {
    public:
        Field(const index_type nMax, const int nDim = 1, const std::string name = "null", const std::string units = "n/a");
        
        inline std::string name() const {return _name;}
        inline std::string units() const {return _units;}
        
        inline index_type nMax() const {return _nMax;}
        inline index_type n() const {return comp0.size();}
        inline int nDim() const {return _nDim;}
        
        inline scalar_type getScalar(const index_type ind) const {return comp0[ind];}
        inline XyzVector get2dVector(const index_type ind) const {return XyzVector(comp0[ind], comp1[ind]);}
        inline XyzVector get3dVector(const index_type ind) const {return XyzVector(comp0[ind], comp1[ind], comp2[ind]);}
        
        void insert(const scalar_type fx);
        void insert(const scalar_type fx, const scalar_type fy, const scalar_type fz = 0.0);
        void insert(const XyzVector& vec);
        
        void replace(const index_type ind, const scalar_type fx);
        void replace(const index_type ind, const scalar_type fx, const scalar_type fy, const scalar_type fz = 0.0);
        void replace(const index_type ind, const XyzVector& vec);
        
        void initializeToConstant(const scalar_type val = 0.0);
        void initializeToScalarFunction(const Coords* crds, const AnalyticFunction* fn);
        void initializeToVectorFunction(const Coords* crds, const AnalyticFunction* fn);
    
    protected:
        std::vector<scalar_type> comp0;
        std::vector<scalar_type> comp1;
        std::vector<scalar_type> comp2;
        
        std::string _name;
        std::string _units;
        
        index_type _nMax;
        int _nDim;
        
        static std::unique_ptr<Logger> log;
};

}

#endif