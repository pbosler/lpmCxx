#ifndef _LPM_FACES_H_
#define _LPM_FACES_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmCoords.h"
#include "LpmEdges.h"
#include "LpmXyzVector.h"
#include "LpmLogger.h"
#include <vector>
#include <tuple>
#include <map>
#include <memory>
#include <string>

namespace Lpm {

class Faces {
    public:
        typedef std::tuple<index_type, index_type, index_type, index_type> quad_index_type;
        virtual ~Faces() {};
        
        XyzVector centroid(const index_type i, const bool lagrangian = false) const;
        
        inline bool isDivded(const index_type i) const {return _hasChildren[i];}
        inline bool hasChildren(const index_type i) const {return _hasChildren[i];}
        inline bool isLeaf(const index_type i) const {return !_hasChildren[i];}
        
        inline index_type positiveCell(const index_type i) const {return _positiveCell[i];}
        inline index_type negativeCell(const index_type i) const {return _negativeCell[i];}
        inline void setPositiveCell(const index_type i, const index_type cellInd) {_positiveCell[i] = cellInd;}
        inline void setNegativeCell(const index_type i, const index_type cellInd) {_negativeCell[i] = cellInd;}
        
        inline bool is3d() const {return _is3d;}
        
        inline index_type n() const {return _edgeInds.size();}
        inline index_type nMax() const {return _nMax;}
        
        index_type nDivided() const;
        inline index_type nLeaves() const {return n() - nDivided();}
        
        inline std::vector<index_type> edgeIndices(const index_type i) const {return _edgeInds[i];}
        std::vector<index_type> vertexIndices(const index_type i) const ;
        
        inline scalar_type area(const index_type i) const {return _area[i];}
        void setArea(const index_type i, const scalar_type nA);
        
        inline quad_index_type children(const index_type i) const {return _children[i];}
        
        void insert(const std::vector<index_type>& edgeInds);
        
        virtual void divide(const index_type i) = 0;
    
        scalar_type computeArea(const index_type i);
        
        scalar_type surfaceArea() const;
        
        void resetAreas();
        
        inline bool edgeIsPositive(const index_type faceInd, const index_type edgeInd) const {
            return (edges->leftFace(edgeInd) == faceInd); }
            
        bool verifyConnectivity(const index_type i) const;
            
    protected:
        Faces(const index_type nMax, const index_type nMaxEdgesPerFace, 
            const std::shared_ptr<Edges> edge_ptr, const std::shared_ptr<Coords> crd_ptr,
            const std::shared_ptr<Coords> lag_crd_ptr = 0,  
            const bool sim3d = false);
        
        std::map<std::string, std::unique_ptr<Field>> fieldMap;
        
        std::vector<std::vector<index_type>> _edgeInds;
        std::vector<scalar_type> _area;
        std::vector<bool> _hasChildren;
        std::vector<index_type> _parent;
        std::vector<quad_index_type> _children;
        std::vector<index_type> _positiveCell;
        std::vector<index_type> _negativeCell;
        
        std::shared_ptr<Edges> edges;
        std::shared_ptr<Coords> crds;
        std::shared_ptr<Coords> lagCrds;
        
        index_type _nMax;
        index_type _nMaxEdges;
        bool _is3d;
        
        static std::unique_ptr<Logger> log;
};

}

#endif