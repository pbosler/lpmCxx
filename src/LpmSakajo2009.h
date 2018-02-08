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
    typedef std::map<MultiIndex, scalar_type> coeff_type;
    typedef std::map<MultiIndex, std::vector<scalar_type>> moment_type;
    
    SakajoNode(const Box3d& bbox, SakajoNode* parent, const int max_series_order);
    
    std::vector<MultiIndex> getKeys() const;
    
    bool multipoleAcceptance(const XyzVector& queryVec, const scalar_type meshSize, const scalar_type nu) const;
    
    std::string coeffString() const;
    std::string momentString() const;
    
    std::string infoString() const override;
    
    void calc_moments(const XyzVector& srcY, const scalar_type srcGamma);
    void calc_coeffs(const XyzVector& tgtVec, const scalar_type sph_rad = 1.0, const scalar_type smoother = 0.0);
    
    coeff_type coeffs;
    moment_type moments;
    
};

class SakajoTree : public Tree {
    public :
        SakajoTree(const int max_series_order, const int max_tree_depth, const scalar_type sphere_radius = 1.0, 
            const scalar_type smooth_param = 0.0, const int prank = 0);
        
        void buildTree(const TREE_DEPTH_CONTROL depth_type, const index_type intParam, const bool do_shrink=false) override {};
        
        void computeMoments(const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
        
        void computeMoments(const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> vorticity, const std::shared_ptr<Field> area);
        
        void writeToVtk(const std::string& filename, const std::string& desc = "") const override;
        
        SakajoNode* getRoot() {return dynamic_cast<SakajoNode*>(_root.get());}
        
        void computeVelocity(std::shared_ptr<Field> outputVelocity, const std::shared_ptr<SphericalCoords> crds, 
                const std::shared_ptr<Field> circ, const scalar_type meshSize, const scalar_type nuPower);
        
//         XyzVector computeVelocity(const index_type tgt_ind, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
        
//         XyzVector computeVelocity(const index_type tgt_ind, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> vorticity, const std::shared_ptr<Field> area);
//         
//         scalar_type computeStreamFunction(const index_type tgt_ind, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
        
        void printAll() const;
        
    protected:
        //bool multipoleAcceptance(const SakajoNode* node, const scalar_type meshSize, const scalar_type nuPower, const XyzVector& queryVec) const;
    
        void generateTree(SakajoNode* node, const int j);
        
        void printNodeInfo(const SakajoNode* node) const;
        
        void nodeMoments(SakajoNode* node, const int k, const XyzVector vecy, const index_type yind, const scalar_type Gamma);
        
        void velocity(XyzVector& vel, SakajoNode* node, const int k, const XyzVector& tgtVec, const index_type tgtInd,
            const scalar_type meshSize, const scalar_type nuPower, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
        
        void streamFn(scalar_type& psi, SakajoNode* node, const int k, const XyzVector& tgtVec, 
            const scalar_type meshSize, const scalar_type nuPower, const std::shared_ptr<SphericalCoords> crds, const std::shared_ptr<Field> circ);
        
        
        
        scalar_type greens(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type smooth_param = 0.0) const;
        
        int _maxSeriesOrder;
        int _maxTreeDepth;
        scalar_type _sphRadius;
        scalar_type _smooth;

};

XyzVector biotSavart(const XyzVector& tgtVec, const XyzVector& srcVec, const scalar_type sphRadius = 1.0, const scalar_type smooth_param = 0.0);

}
#endif