

#add_library( lpmLib GlobalConstants.cpp OutputMessage.cpp Logger.cpp
#					LpmXyzVector.hpp LpmCoords.hpp
#					LpmParticles.hpp LpmScalarField.hpp LpmField2d.hpp LpmField3d.hpp
#				    #Edges.cpp Faces.cpp PolyMesh2d.cpp 
#				    #RefinementType.cpp 
#				    )
				    
add_library( lpm LpmXyzVector.cpp LpmOutputMessage.cpp LpmLogger.cpp LpmCoords.cpp 
             LpmEuclideanCoords.cpp LpmSphericalCoords.cpp LpmField.cpp LpmEdges.cpp
             LpmFaces.cpp LpmTriFaces.cpp LpmQuadFaces.cpp LpmMeshSeed.cpp
)				    
				    
#add_library( oldLpmLib GlobalConstants.cpp OutputMessage.cpp Logger.cpp xyzVector.cpp
#					   Particles.cpp Edges.cpp Faces.cpp PolyMesh2d.cpp Field.cpp )				    
include_directories( "${PROJECT_SRC_DIR}/src" )				    
				 
#
# install targets
#				 
#install ( FILES GlobalConstants.h OutputMessage.h Logger.h 
#				LpmXyzVector.hpp LpmCoords.hpp LpmSphereCoords.hpp
#				LpmParticles.hpp LpmScalarField.hpp LpmField2d.hpp LpmField3d.hpp
#				#Edges.h Faces.h PolyMesh2d.h
#				#RefinementType.h
#		  DESTINATION include )

install (FILES LpmXyzVector.h LpmOutputMessage.h LpmLogger.h LpmCoords.h LpmEuclideanCoords.h LpmSphericalCoords.h
               LpmField.h LpmAnalyticFunctions.h LpmEdges.h LpmFaces.h LpmTriFaces.h LpmQuadFaces.h LpmMeshSeed.h
         DESTINATION include)
		  
install ( TARGETS lpm DESTINATION lib )  
		  
