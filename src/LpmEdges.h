#ifndef _LPM_EDGES_H_
#define _LPM_EDGES_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmCoords.h"
#include "LpmXyzVector.h"
#include "LpmLogger.h"
#include "LpmOutputMessage.h"
#include "LpmField.h"
#include <vector>
#include <tuple>
#include <memory>
#include <map>

namespace Lpm {

class Edges {
    public:
        typedef std::pair<index_type, index_type> index_pair_type;
        
        Edges(const index_type nMax, const std::shared_ptr<Coords> crd_ptr, const std::shared_ptr<Coords> lag_crd_ptr = 0);
        virtual ~Edges() {};
        
        inline void makeLagrangian(const std::shared_ptr<Coords> lag_crd_ptr) {lagCrds = lag_crd_ptr;}
    
        inline index_type nMax() const {return _nMax;}
        inline index_type n() const {return _orig.size();}
        inline index_type size() const {return _orig.size();}
        inline index_type orig(const index_type i) const {return _orig[i];}
        inline index_type dest(const index_type i) const {return _dest[i];}
        inline index_type rightFace(const index_type i) const {return _rightFace[i];}
        inline index_type leftFace(const index_type i) const {return _leftFace[i];}
        
        inline index_pair_type children(const index_type i) const {return index_pair_type(_child0[i], _child1[i]);}
        inline index_pair_type vertices(const index_type i) const {return index_pair_type(_orig[i], _dest[i]);}
        inline bool isDivided(const index_type i) const {return _hasChildren[i];}
        inline bool hasChildren(const index_type i) const {return _hasChildren[i];}
        inline bool isLeaf(const index_type i) const {return !_hasChildren[i];}
        
        index_type nDivided() const; 
        index_type nLeaves() const; 
        
        inline bool onBoundary(const index_type i) const {return (_leftFace[i] < 0 || _rightFace[i] < 0);}
        
//         std::vector<index_type> getLeafEdgesFromParent(const index_type i) const;
//         std::vector<index_type> getAllVertices(const index_type i) const;
        
        scalar_type length(const index_type i, const bool lagrangian = false) const;
        scalar_type minLength(const bool lagrangian = false) const;
        scalar_type maxLength(const bool lagrangian = false) const;
        XyzVector midpoint(const index_type i, const bool lagrangian = false) const;
        XyzVector origCoord(const index_type i, const bool lagrangian = false) const;
        XyzVector destCoord(const index_type i, const bool lagrangian = false) const;
        XyzVector edgeVector(const index_type i, const bool lagrangian = false) const;
        
        inline void setLogProc(const int rank) {log->setProcRank(rank);}
        
        inline void setLeftFace(const index_type i, const index_type faceInd) {_leftFace[i] = faceInd;}
        inline void setRightFace(const index_type i, const index_type faceInd) {_rightFace[i] = faceInd;}
        
        virtual void divide(const index_type i);
        virtual void insert(const index_type origInd, const index_type destInd, 
            const index_type leftInd, const index_type rightInd);
            
        std::string edgeRecord(const index_type i) const;
    protected:
        std::vector<index_type> _orig;
        std::vector<index_type> _dest;
        std::vector<index_type> _leftFace;
        std::vector<index_type> _rightFace;
        std::vector<index_type> _child0;
        std::vector<index_type> _child1;
        std::vector<bool> _hasChildren;
        std::vector<index_type> _parent;
        
        index_type _nMax;
        
        std::weak_ptr<Coords> crds;
        std::weak_ptr<Coords> lagCrds;
        
        static std::unique_ptr<Logger> log;
};

}
#endif