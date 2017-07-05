include(CTest)

#add_executable ( lpmBaseObsTest.exe baseTest.cpp )
#target_link_libraries( lpmBaseObsTest.exe lpmLib )
#
#add_executable ( lpmCoordsUnitTest.exe lpmCoordsUnitTest.cpp )
#target_link_libraries( lpmCoordsUnitTest.exe lpmLib )
#
#add_executable ( lpmParticlesAndFieldsUnitTest.exe lpmParticlesAndFieldTest.cpp )
#target_link_libraries( lpmParticlesAndFieldsUnitTest.exe lpmLib )

#add_executable ( lpmParticlesAndFieldTest.exe ParticlesAndFieldTest.cpp )
#target_link_libraries( lpmParticlesAndFieldTest.exe lpmLib )

#add_executable( lpmPlanarMeshTest.exe PlanarMeshTest.cpp )
#target_link_libraries( lpmPlanarMeshTest.exe lpmLib )

#add_executable( lpmSphereMeshTest.exe sphereVTKMeshTest.cpp)
#target_link_libraries( lpmSphereMeshTest.exe oldLpmLib )

#install ( TARGETS lpmBaseObsTest.exe lpmParticlesAndFieldTest.exe lpmPlanarMeshTest.exe lpmSphereMeshTest.exe
#		  DESTINATION bin )

ADD_EXECUTABLE(lpmXyzVectorUnitTest.exe LpmXyzVectorUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmXyzVectorUnitTest.exe lpm)
add_test( xyzVecTest lpmXyzVectorUnitTest.exe)

ADD_EXECUTABLE(lpmOutputMessageUnitTest.exe LpmOutputMessageUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmOutputMessageUnitTest.exe lpm)
add_test( outputMessageTest lpmOutputMessageUnitTest.exe)

ADD_EXECUTABLE(lpmLoggerUnitTest.exe LpmLoggerUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmLoggerUnitTest.exe lpm)
add_test( loggerTest COMMAND mpirun -np 2 lpmLoggerUnitTest.exe)

ADD_EXECUTABLE(lpmCoordsUnitTest.exe LpmCoordsUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmCoordsUnitTest.exe lpm)
add_test( coordsTest lpmCoordsUnitTest.exe)

ADD_EXECUTABLE(lpmFieldUnitTest.exe LpmFieldUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmFieldUnitTest.exe lpm)
add_test( fieldTest lpmFieldUnitTest.exe)

ADD_EXECUTABLE(lpmEdgesUnitTest.exe LpmEdgesUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmEdgesUnitTest.exe lpm)
add_test( edgeTest lpmEdgesUnitTest.exe)

ADD_EXECUTABLE(lpmFacesUnitTest.exe LpmFacesUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmFacesUnitTest.exe lpm)
add_test( faceTest lpmFacesUnitTest.exe)

ADD_EXECUTABLE(lpmMeshSeedUnitTest.exe LpmMeshSeedUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmMeshSeedUnitTest.exe lpm)
add_test( meshSeedUnitTest lpmMeshSeedUnitTest.exe)

ADD_EXECUTABLE(lpmPolyMeshUnitTest.exe LpmPolyMesh2dUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmPolyMeshUnitTest.exe lpm)
ADD_TEST(polymesh2dUnitTest lpmPolyMeshUnitTest.exe)

ADD_EXECUTABLE(lpmMPIReplicatedDataUnitTest.exe LpmMPIReplicatedDataUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmMPIReplicatedDataUnitTest.exe lpm)
ADD_TEST(mpiReplicatedDataUnitTest COMMAND mpirun -np 2 lpmMPIReplicatedDataUnitTest.exe)

ADD_EXECUTABLE(lpmTimerUnitTest.exe LpmTimerUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmTimerUnitTest.exe lpm)
ADD_TEST(timerUnitTest COMMAND mpirun -np 2 lpmTimerUnitTest.exe)

ADD_EXECUTABLE(lpmPoissonSolver.exe LpmPoissonSolverTests.cpp)
TARGET_LINK_LIBRARIES(lpmPoissonSolver.exe lpm)
ADD_TEST(poissonSolver lpmPoissonSolver.exe)
SET_TESTS_PROPERTIES(poissonSolver PROPERTIES TIMEOUT 300)

#add_test( baselibTest lpmBaseObsTest.exe )
#add_test( coordsTest lpmCoordsUnitTest.exe )
#add_test( particlesAndFieldTest lpmParticlesAndFieldsUnitTest.exe )
#add_test( particlesAndFieldTest lpmParticlesAndFieldTest.exe )
#add_test( planarMeshTest lpmPlanarMeshTest.exe )
#add_test( sphereMeshTest lpmSphereMeshTest.exe )
