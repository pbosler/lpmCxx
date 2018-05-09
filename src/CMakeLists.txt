#
#   build
#	    
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src" "${PROJECT_BINARY_DIR}")

add_library( lpm LpmXyzVector.cpp LpmOutputMessage.cpp LpmLogger.cpp LpmCoords.cpp 
             LpmEuclideanCoords.cpp LpmSphericalCoords.cpp LpmField.cpp LpmEdges.cpp
             LpmFaces.cpp LpmTriFaces.cpp LpmQuadFaces.cpp LpmMeshSeed.cpp LpmPolyMesh2d.cpp
             LpmVtkFileIO.cpp LpmParticles.cpp LpmMeshedParticles.cpp LpmMPIReplicatedData.cpp
             LpmTimer.cpp LpmDirectSum.cpp
             LpmBox3d.cpp LpmOctree.cpp LpmSakajo2009.cpp
             LpmTreeSum.cpp 
)


#
#   install 
#
install (FILES LpmConfig.h LpmXyzVector.h LpmOutputMessage.h LpmLogger.h LpmCoords.h LpmEuclideanCoords.h LpmSphericalCoords.h
               LpmField.h LpmAnalyticFunctions.h LpmEdges.h LpmFaces.h LpmTriFaces.h LpmQuadFaces.h LpmMeshSeed.h
               LpmPolyMesh2d.h LpmVtkFileIO.h LpmParticles.h LpmMeshedParticles.h LpmMPIReplicatedData.h
               LpmTimer.h LpmMultiIndex.h LpmSakajo2009.h LpmTreeSum.h 
               LpmDirectSum.h 
               LpmAosTypes.hpp
         DESTINATION include)
		  
install ( TARGETS lpm DESTINATION lib )  
		  
