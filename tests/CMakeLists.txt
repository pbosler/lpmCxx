include(CTest)

add_executable ( lpmBaseObsTest.exe baseTest.cpp )
target_link_libraries( lpmBaseObsTest.exe lpmLib )

add_executable ( lpmCoordsUnitTest.exe lpmCoordsUnitTest.cpp )
target_link_libraries( lpmCoordsUnitTest.exe lpmLib )

add_executable ( lpmParticlesAndFieldsUnitTest.exe lpmParticlesAndFieldTest.cpp )
target_link_libraries( lpmParticlesAndFieldsUnitTest.exe lpmLib )

#add_executable ( lpmParticlesAndFieldTest.exe ParticlesAndFieldTest.cpp )
#target_link_libraries( lpmParticlesAndFieldTest.exe lpmLib )

#add_executable( lpmPlanarMeshTest.exe PlanarMeshTest.cpp )
#target_link_libraries( lpmPlanarMeshTest.exe lpmLib )

#add_executable( lpmSphereMeshTest.exe sphereVTKMeshTest.cpp)
#target_link_libraries( lpmSphereMeshTest.exe lpmLib )

#install ( TARGETS lpmBaseObsTest.exe lpmParticlesAndFieldTest.exe lpmPlanarMeshTest.exe lpmSphereMeshTest.exe
#		  DESTINATION bin )

add_test( baselibTest lpmBaseObsTest.exe )
add_test( coordsTest lpmCoordsUnitTest.exe )
add_test( particlesAndFieldTest lpmParticlesAndFieldsUnitTest.exe )
#add_test( particlesAndFieldTest lpmParticlesAndFieldTest.exe )
#add_test( planarMeshTest lpmPlanarMeshTest.exe )
#add_test( sphereMeshTest lpmSphereMeshTest.exe )
