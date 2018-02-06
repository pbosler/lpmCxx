#ifndef _LPM_SAKAJO_2009_H_
#define _LPM_SAKAJO_2009_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmXyzVector.h"
#include "LpmSphericalCoords.h"
#include "LpmField.h"
#include "LpmOctree.h"
#include "LpmLogger.h"
#include "LpmBox3d.h"
#include "LpmMultiIndex.h"
#include "LpmXyzVector.h"
#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <map>

namespace Lpm {

struct SakajoNode : public Node {
    typedef std::map<MultiIndex, std::vector<scalar_type>> coeff_type;
    
    SakajoNode(const Box3d& bbox, SakajoNode* parent, const int max_series_order);
    
    coeff_type coeffs;
    coeff_type moments;
};

class SakajoTree : public Tree {
    public :
        SakajoTree(const int max_series_order, const int max_tree_depth, const scalar_type sphere_radius = 1.0, 
            const scalar_type smooth_param = 0.0, const int prank = 0);
        
        void buildTree(const TREE_DEPTH_CONTROL depth_type, const index_type intParam, const bool do_shrink=false) override {};
        
        void computeCoefficients(const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
        
        void computeCoefficients(const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> vorticity, const std::shared_ptr<Field> area);
        
//         XyzVector computeVelocity(const index_type tgt_ind, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
//         XyzVector computeVelocity(const index_type tgt_ind, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> vorticity, const std::shared_ptr<Field> area);
//         
//         scalar_type computeStreamFunction(const index_type tgt_ind, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
        
    protected:
        bool multipoleAcceptance(const SakajoNode* node, const scalar_type meshSize, const scalar_type nuPower, const XyzVector& queryVec) const;
    
        void generateTree(SakajoNode* node, const int j);
        
        void nodeCoeffs(SakajoNode* node, const XyzVector& tgtVec);
        void nodeMoments(SakajoNode* node, const int k, const XyzVector vecy, const index_type yind, const scalar_type Gamma);
        
        void velocity(XyzVector& vel, SakajoNode* node, const int k, const XyzVector& tgtVec, 
            const scalar_type meshSize, const scalar_type nuPower, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
        
        void streamFn(scalar_type& psi, SakajoNode* node, const int k, const XyzVector& tgtVec, 
            const scalar_type meshSize, const scalar_type nuPower, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
        
        XyzVector biotSavart(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type smooth_param = 0.0) const;
        
        scalar_type greens(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type smooth_param = 0.0) const;
        
        int _maxSeriesOrder;
        int _maxTreeDepth;
        scalar_type _sphRadius;
        scalar_type _smooth;

};

}
#endif