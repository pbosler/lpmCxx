#ifndef _LPM_COORDS_BASE_H_
#define _LPM_COORDS_BASE_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmLogger.h"
#include <memory>
#include <vector>

namespace Lpm {

class Coords {
    public:
        virtual ~Coords() {};
        virtual scalar_type distance(const index_type indA, const index_type indB) const = 0;
        virtual scalar_type distance(const XyzVector& v0, const index_type ind) const = 0;
        virtual scalar_type distance(const XyzVector& v0, const XyzVector& v1) const = 0;
        virtual XyzVector midpoint(const index_type indA, const index_type indB) const = 0;
        virtual XyzVector centroid(const std::vector<index_type>& inds) const = 0;
        virtual scalar_type triArea(const index_type indA, const index_type indB, const index_type indC) const = 0;
        virtual scalar_type triArea(const XyzVector& v0, const index_type indA, const index_type indB) const = 0;
        
        XyzVector crossProduct(const index_type indA, const index_type indB) const;
        scalar_type dotProduct(const index_type indA, const index_type indB) const;
        
        std::vector<XyzVector> getVectors(const std::vector<index_type> inds) const;
        
        inline index_type nMax() const {return _nMax;}
        inline index_type n() const {return x.size();}
        
        inline GeometryType geometry() const {return _geometry;}
        inline bool is2d() const {return (_geometry == PLANAR_GEOMETRY);}
        
        void scaleAll(const scalar_type multiplier);
        void normalizeAll();
        scalar_type magnitude(const index_type ind) const;
        
        void replace(const index_type ind, const scalar_type nx, const scalar_type ny, const scalar_type nz = 0.0);
        void replace(const index_type ind, const XyzVector& vec);
        
        void insert(const scalar_type nx, const scalar_type ny, const scalar_type nz = 0.0);
        void insert(const XyzVector& vec);
        
        inline XyzVector getVec(const index_type ind) const {return XyzVector(x[ind], y[ind], (_geometry == PLANAR_GEOMETRY ? 0.0 : z[ind]));}
        
        std::string listAllCoords() const;
        
    protected:
        Coords(const index_type nMax, const bool coords3d = true);
        
        std::vector<scalar_type> x;
        std::vector<scalar_type> y;
        std::vector<scalar_type> z;
        int _nMax;
        GeometryType _geometry;
        
        static std::unique_ptr<Logger> log;
};

}

std::ostream& operator << (std::ostream& os, const Lpm::Coords& crds);

#endif