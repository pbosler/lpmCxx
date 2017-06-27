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
        
        XyzVector centroid(const index_type i) const;
        
        inline bool isDivded(const index_type i) const {return _hasChildren[i];}
        inline bool hasChildren(const index_type i) const {return _hasChildren[i];}
        inline bool isLeaf(const index_type i) const {return !_hasChildren[i];}
        
        inline std::vector<index_type> edgeIndices(const index_type i) const {return _edgeInds[i];}
        std::vector<index_type> vertexIndices(const index_type i) const ;
        
        inline scalar_type area(const index_type i) const {return _area[i];}
        void setArea(const index_type i, const scalar_type nA);
        
        inline quad_index_type children(const index_type i) const {return _children[i];}
        
        virtual void insert(const std::vector<index_type>& edgeInds);
        
        virtual void divide(const index_type i) = 0;
    
        void computeArea(const index_type i);
        
        inline bool edgeIsPositive(const index_type faceInd, const index_type edgeInd) const {
            return (edges->leftFace(edgeInd) == faceInd); }
            
        bool verifyConnectivity(const index_type i) const;
            
    protected:
        Faces(const index_type nMax, const index_type nMaxEdgesPerFace, 
            const std::shared_ptr<Edges> edge_ptr, const std::shared_ptr<Coords> crd_ptr, 
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
        
        index_type _nMax;
        index_type _nLeaves;
        index_type _nMaxEdges;
        
        static std::unique_ptr<Logger> log;
};

}

#endif