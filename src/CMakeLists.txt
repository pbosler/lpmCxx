#
#   build
#
include_directories( "${PROJECT_SRC_DIR}/src" )				    

add_library( lpm LpmXyzVector.cpp LpmOutputMessage.cpp LpmLogger.cpp LpmCoords.cpp 
             LpmEuclideanCoords.cpp LpmSphericalCoords.cpp LpmField.cpp LpmEdges.cpp
             LpmFaces.cpp LpmTriFaces.cpp LpmQuadFaces.cpp LpmMeshSeed.cpp LpmPolyMesh2d.cpp
             LpmVtkFileIO.cpp
)


#
#   install 
#
install (FILES LpmXyzVector.h LpmOutputMessage.h LpmLogger.h LpmCoords.h LpmEuclideanCoords.h LpmSphericalCoords.h
               LpmField.h LpmAnalyticFunctions.h LpmEdges.h LpmFaces.h LpmTriFaces.h LpmQuadFaces.h LpmMeshSeed.h
               LpmPolyMesh2d.h
         DESTINATION include)
		  
install ( TARGETS lpm DESTINATION lib )  
		  
