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
        inline void rename(const std::string& nname) {_name = nname;}
        
        inline index_type nMax() const {return _nMax;}
        inline index_type n() const {return comp0.size();}
        inline int nDim() const {return _nDim;}
        
        inline scalar_type getScalar(const index_type ind) const {return comp0[ind];}
        inline XyzVector get2dVector(const index_type ind) const {return XyzVector(comp0[ind], comp1[ind]);}
        inline XyzVector get3dVector(const index_type ind) const {return XyzVector(comp0[ind], comp1[ind], comp2[ind]);}
        
        inline void clear() {comp0.clear(); comp1.clear(); comp2.clear();}
        
        inline scalar_type* getPtrToData(const int comp = 0) {
            scalar_type* result = NULL;
            switch (comp) {
                case (0) : {
                    result = &comp0[0];
                    break;
                }
                case (1) : {
                    result = &comp1[0];
                    break;
                }
                case (2) : {
                    result = &comp2[0];
                    break;
                }
            }
            return result;
        }
        
        void insert(const scalar_type fx);
        void insert(const scalar_type fx, const scalar_type fy, const scalar_type fz = 0.0);
        void insert(const XyzVector& vec);
        
        void replace(const index_type ind, const scalar_type fx);
        void replace(const index_type ind, const scalar_type fx, const scalar_type fy, const scalar_type fz = 0.0);
        void replace(const index_type ind, const XyzVector& vec);
        
        void initializeToConstant(const Coords* crds, const scalar_type val = 0.0);
        void initializeToConstant(const index_type nn, const scalar_type val = 0.0);
        void initializeToScalarFunction(const Coords* crds, const AnalyticFunction* fn);
        void initializeToVectorFunction(const Coords* crds, const AnalyticFunction* fn);
        
        void abs();
        void update(const scalar_type a, const std::shared_ptr<Field>& other, 
                    const scalar_type b = 0.0, const std::shared_ptr<Field>& other1 = NULL);
                    
        void scale(const scalar_type multiplier);
        scalar_type maxScalarVal() const;
        scalar_type minScalarVal() const;
        
        inline void setLogProc(const int rank) {log->setProcRank(rank);}
        
        std::string listFieldValues() const;
    
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