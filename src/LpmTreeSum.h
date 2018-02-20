#ifndef _LPM_TREESUM_H_
#define _LPM_TREESUM_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmOctree.h"
#include "LpmMultiIndex.h"
#include "LpmXyzVector.h"
#include "LpmField.h"
#include "LpmCoords.h"
#include "LpmBox3d.h"
#include "LpmScalarKernel.h"
#include "LpmVectorKernel.h"
#include <map>
#include <string>
#include <memory>
#include <vector>


namespace Lpm {

struct SumNode : public Node {
    typedef std::map<MultiIndex, scalar_type> coeff_type;
    typedef std::map<MultiIndex, std::vector<scalar_type>> moment_type;
    
    SumNode(const Box3d& bbox, SumNode* parent, const std::vector<index_type>& crd_indices, const int max_series_order);
    
    std::vector<MultiIndex> getKeys() const;
    
    virtual bool multipoleAcceptance(const XyzVector& tgtVec, const scalar_type tol) const;
    
    std::string coeffString() const;
    std::string momentString() const;
    
    std::string infoString() const;
    
    void calc_moments(const XyzVector& srcVec, const scalar_type srcStrength);
    void calc_coeffs(const XyzVector& tgtVec, const scalar_type smoother = 0.0);
    
    coeff_type coeffs;
    moment_type moments;
};

class TreeSum : public Tree {
    public:
        TreeSum(const std::shared_ptr<Coords> crds, const int max_series_order, const int max_particles_per_box, 
    const scalar_type smoother, const bool do_shrink=false);
        
        inline void initVectorKernel(VectorKernel* kernel){vKern = kernel;}
        inline void initScalarKernel(ScalarKernel* kernel){sKern = kernel;}
        
        void computeMoments(const std::shared_ptr<Field> strength);
        
        void computeScalarSum(std::shared_ptr<Field> outputScalar, const std::shared_ptr<Coords> crds,
            const std::shared_ptr<Field> strength);
        
        void computeVectorSum(std::shared_ptr<Field> outputVector, const std::shared_ptr<Coords> crds,
            const std::shared_ptr<Field> strength);
            
    protected:
        void generateTree(SumNode* node, const Coords* crds);
        
        void nodeMoments(SumNode* node, const XyzVector& srcVec, const scalar_type srcStrength);    
        
        void vecSum(XyzVector& vel, const XyzVector& tgtVec, SumNode* node, const Coords* crds, const scalar_type mp_tol, 
             const std::shared_ptr<Field> strength);
        
        int _maxSeriesOrder;
        int _maxParticlesPerNode;
        scalar_type _smooth;
        bool _shrink;
        
        ScalarKernel* sKern;
        VectorKernel* vKern;
};

}

#endif
