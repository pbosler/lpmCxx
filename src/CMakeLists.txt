add_library( lpmLib GlobalConstants.cpp OutputMessage.cpp Logger.cpp
					LpmXyzVector.hpp LpmCoords.hpp
					LpmParticles.hpp LpmScalarField.hpp LpmField2d.hpp LpmField3d.hpp
				    #Edges.cpp Faces.cpp PolyMesh2d.cpp 
				    #RefinementType.cpp 
				    )
include_directories( "${PROJECT_SRC_DIR}/src" )				    
				 
#
# install targets
#				 
install ( FILES GlobalConstants.h OutputMessage.h Logger.h 
				LpmXyzVector.hpp LpmCoords.hpp LpmSphereCoords.hpp
				LpmParticles.hpp LpmScalarField.hpp LpmField2d.hpp LpmField3d.hpp
				#Edges.h Faces.h PolyMesh2d.h
				#RefinementType.h
		  DESTINATION include )
		  
install ( TARGETS lpmLib DESTINATION lib )  
		  
