#ifndef _LPM_TREE_SUM_H_
#define _LPM_TREE_SUM_H_

#include "LpmConfig.h"
#include "LpmTypeDefs.h"
#include "LpmLogger.h"
#include "LpmOctree.h"
#include "LpmXyzVector.h"
#include "LpmCoords.h"
#include "LpmField.h"
#include "LpmMultiIndex.h"
#include "LpmTaylorSeries3d.h"
#include "LpmScalarKernel.h"
#include "LpmMeshedParticles.h"
#include "LpmMPIReplicatedData.h"
#include <vector>
#include <string>
#include <cmath>
#include <memory>

namespace Lpm {

struct SumNode : public Node {
    typedef std::unique_ptr<TaylorSeries3d> series_ptr_type;

    SumNode(const Box3d& bbox, Node* pparent = NULL, const std::vector<index_type>& crdInds = std::vector<index_type>(),
         const int maxSeriesOrder = 0, ScalarKernel* kernel = NULL);

    series_ptr_type series;
    
    void computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcStrength);
    void computeMoments(const std::shared_ptr<Coords> crds, const std::shared_ptr<Field> srcVals, 
        const std::shared_ptr<Field> srcWeights);
    void computeMoments(const std::shared_ptr<Faces> faces, const std::shared_ptr<Field> srcVals);
    void computeCoeffs(const XyzVector& tgtVec, const scalar_type param = 0.0);
    
    scalar_type seriesSum() const {return series->sum();}
    
    inline bool isFar(const XyzVector& tgtVec, const int maxTreeDepth, const scalar_type nuParam = 1.0) const {
        const scalar_type hnu = std::pow(std::pow(2.0, -maxTreeDepth), nuParam);
        const XyzVector vec = tgtVec - box.centroid();
        return box.radius() <= hnu * vec.magnitude();
    }
    
    inline bool isNear(const XyzVector& tgtVec, const int maxTreeDepth, const scalar_type nuParam) const {
        return !isFar(tgtVec, maxTreeDepth, nuParam);
    }
    
    bool momentsReady;
};

class TreeSum : public Tree {
    public:
        TreeSum(const std::shared_ptr<Coords> crds, const scalar_type maxAspectRatio, 
            const std::shared_ptr<ScalarKernel> kernel, const int maxSeriesOrder, const scalar_type seriesParam = 0.0,
            const int prank = 0, const scalar_type farFieldParam = 1.0);
        
        TreeSum(std::shared_ptr<MeshedParticles> pmesh, const std::shared_ptr<ScalarKernel> kernel, 
            const std::string tgtFieldName, const std::string srcFieldName, const scalar_type maxAspectRatio, 
            const int maxSeriesOrder, const scalar_type seriesParam, const int prank = 0, 
            const scalar_type farFieldParam = 1.0);

        inline void resetMoments() {setRecomputeMomentsTrue(_root.get());}

        void buildTree(const index_type maxCoordsPerNode) override;
        
        void buildTreeFromMesh(const index_type maxCoordsPerNode);
        
        scalar_type computeSum(const XyzVector& tgtLoc, const std::shared_ptr<Field> srcStrength);
        scalar_type computeSum(const XyzVector& tgtLoc, const std::shared_ptr<Field> srcVals, 
            const std::shared_ptr<Field> srcWeights);

        void meshSolve(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces);
        void meshBroadcast(const MPIReplicatedData& mpiVerts, const MPIReplicatedData& mpiFaces) const;

    protected:
        std::weak_ptr<ScalarKernel> _kernel;
        scalar_type _nuParam;
        int _maxP;
        scalar_type _seriesParam;
        
        void recursiveSum(scalar_type& sum, const XyzVector& tgtLoc, const std::shared_ptr<Field> srcVals, 
            const std::shared_ptr<Field> srcWeights, SumNode* node, const index_type& tree_depth);
            
        void recursiveSum(scalar_type& sum, const index_type& tgtInd, const std::shared_ptr<Field> srcVals, 
            const std::shared_ptr<Field> srcWeights, SumNode* node, const index_type& tree_depth);
        
        scalar_type recursiveSumMeshVertices(const index_type& tgtInd, std::shared_ptr<Field> srcVals, 
            SumNode* node, const index_type& tree_depth); 
        
        scalar_type recursiveSumMeshFaces(const index_type& tgtInd, std::shared_ptr<Faces> faces, std::shared_ptr<Field> srcVals,
            SumNode* node, const index_type& tree_depth);
        
        void setRecomputeMomentsTrue(Node* node);
    
        void generateTree(Node* node, const index_type maxCoordsPerNode) override;
        
        void generateTreeFromMesh(Node* node, const index_type maxCoordsPerNode);
        
        // mesh solver variables
        std::weak_ptr<PolyMesh2d> _mesh;
        std::weak_ptr<Field> _vertTgt;
        std::weak_ptr<Field> _faceTgt;
        std::weak_ptr<Field> _faceSrc;
        
        // meshfree solver variables
        std::weak_ptr<Coords> _srcLocs;
        std::weak_ptr<Field> _srcVals;
        std::weak_ptr<Field> _tgtVals;
        std::weak_ptr<Field> _srcWeights;
};


}

#endif
