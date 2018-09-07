include(CTest)

include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src" "${PROJECT_BINARY_DIR}" ${VTK_INCLUDE_DIRS})

ADD_EXECUTABLE(lpmXyzVectorUnitTest.exe LpmXyzVectorUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmXyzVectorUnitTest.exe lpm)
add_test( xyzVecTest lpmXyzVectorUnitTest.exe)

ADD_EXECUTABLE(lpmOutputMessageUnitTest.exe LpmOutputMessageUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmOutputMessageUnitTest.exe lpm)
add_test( outputMessageTest lpmOutputMessageUnitTest.exe)

ADD_EXECUTABLE(lpmLoggerUnitTest.exe LpmLoggerUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmLoggerUnitTest.exe lpm)
add_test(NAME loggerTest COMMAND mpirun -np 2 lpmLoggerUnitTest.exe)

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
ADD_TEST(NAME mpiReplicatedDataUnitTest COMMAND mpirun -np 2 lpmMPIReplicatedDataUnitTest.exe)

ADD_EXECUTABLE(lpmTimerUnitTest.exe LpmTimerUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmTimerUnitTest.exe lpm)
ADD_TEST(NAME timerUnitTest COMMAND mpirun -np 2 lpmTimerUnitTest.exe)

ADD_EXECUTABLE(lpmPoissonSolver.exe LpmPoissonSolverTests.cpp)
TARGET_LINK_LIBRARIES(lpmPoissonSolver.exe lpm)
ADD_TEST(NAME poissonSolver COMMAND mpirun -np 4 lpmPoissonSolver.exe)
SET_TESTS_PROPERTIES(poissonSolver PROPERTIES TIMEOUT 300)

ADD_EXECUTABLE(lpmBox3dUnitTests.exe LpmBox3dUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmBox3dUnitTests.exe lpm)
ADD_TEST(NAME box3dUnitTest COMMAND mpirun -np 1 lpmBox3dUnitTests.exe)

ADD_EXECUTABLE(lpmOctreeUnitTest.exe LpmOctreeUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmOctreeUnitTest.exe lpm)
ADD_TEST(octreeUnitTest lpmOctreeUnitTest.exe)

#ADD_EXECUTABLE(lpmTaylorSeries3dUnitTest.exe LpmTaylorSeries3dUnitTests.cpp)
#TARGET_LINK_LIBRARIES(lpmTaylorSeries3dUnitTest.exe lpm)
#ADD_TEST(NAME taylorSeries3dTest COMMAND mpirun -np 1 lpmTaylorSeries3dUnitTest.exe)

#ADD_EXECUTABLE(lpmTreeSumTest.exe LpmTreeSumTests.cpp)
#TARGET_LINK_LIBRARIES(lpmTreeSumTest.exe lpm)
#ADD_TEST(NAME treeSumTest COMMAND mpirun -np 1 lpmTreeSumTest.exe)

ADD_EXECUTABLE(lpmBiotSavartTest.exe LpmBiotSavartRHWaveTests.cpp)
TARGET_LINK_LIBRARIES(lpmBiotSavartTest.exe lpm)
ADD_TEST(NAME biotSavartTest COMMAND mpirun -np 4 lpmBiotSavartTest.exe)

ADD_EXECUTABLE(lpmMultiIndexUnitTest.exe LpmMultiIndexUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmMultiIndexUnitTest.exe lpm)
ADD_TEST(NAME multiIndexTest COMMAND mpirun -np 1 lpmMultiIndexUnitTest.exe)

ADD_EXECUTABLE(lpmSakajoUnitTest.exe LpmSakajoUnitTests.cpp)
TARGET_LINK_LIBRARIES(lpmSakajoUnitTest.exe lpm)
ADD_TEST(NAME sakajoUnitTest COMMAND mpirun -np 1 lpmSakajoUnitTest.exe)

ADD_EXECUTABLE(lpmAosTypesUnitTest.exe LpmAosUnitTest.cpp)
TARGET_LINK_LIBRARIES(lpmAosTypesUnitTest.exe lpm)
add_test(aosTypeUnitTest lpmAosTypesUnitTest.exe)

ADD_EXECUTABLE(lpmAosParticleUnitTest.exe LpmAosParticleUnitTest.cpp)
TARGET_LINK_LIBRARIES(lpmAosParticleUnitTest.exe lpm)
add_test(aosParticleUnitTest lpmAosParticleUnitTest.exe)

ADD_EXECUTABLE(lpmAosEdgeUnitTest.exe LpmAosEdgeUnitTest.cpp)
TARGET_LINK_LIBRARIES(lpmAosEdgeUnitTest.exe lpm)
add_test(aosEdgeUnitTest lpmAosEdgeUnitTest.exe)

ADD_EXECUTABLE(lpmAosFaceUnitTest.exe LpmAosFaceUnitTest.cpp)
TARGET_LINK_LIBRARIES(lpmAosFaceUnitTest.exe lpm)
add_test(aosFaceUnitTest lpmAosFaceUnitTest.exe)

#ADD_EXECUTABLE(lpmAosMeshSeedTest.exe LpmAosMeshSeedUnitTest.cpp)
#TARGET_LINK_LIBRARIES(lpmAosMeshSeedTest.exe lpm)
#add_test(aosMeshSeedTest lpmAosMeshSeedTest.exe)

#ADD_EXECUTABLE(lpmAosPolyMesh2dUnitTest.exe LpmAosPolyMeshUnitTest.cpp)
#TARGET_LINK_LIBRARIES(lpmAosPolyMesh2dUnitTest.exe lpm)
#add_test(aosPolyMeshTest lpmAosPolyMesh2dUnitTest.exe)


