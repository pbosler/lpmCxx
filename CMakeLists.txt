cmake_minimum_required ( VERSION 2.8 )

project ( LPM )
set (LPM_VERSION_MAJOR 2)
set (LPM_VERSION_MINOR 0)

if ("{$CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
set( CMAKE_CXX_FLAGS_DEBUG, "-g -fopenmp")
set( CMAKE_CXX_FLAGS_RELEASE, "-O2 -fopenmp -Wl,-stack_size -Wl,1000000")
elseif ("{$CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
set( CMAKE_CXX_FLAGS_DEBUG, "-g -fopenmp")
set( CMAKE_CXX_FLAGS_RELEASE, "-O2 -fopenmp")
elseif ("{$CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
set( CMAKE_CXX_FLAGS_DEBUG, "-g -openmp")
set( CMAKE_CXX_FLAGS_RELEASE, "-O2 -openmp")
endif() 

#
#   CMake Checks
#
IF (CMAKE_VERSION VERSION_LESS "3.1")
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    SET (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
  ENDIF ()
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    SET (CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")
  ENDIF ()
ELSE ()
  SET (CMAKE_CXX_STANDARD 11) 
ENDIF ()

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SRC_DIR}/src" "${MPI_INCLUDE_PATH}" "${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src" "${MPI_INCLUDE_PATH}" "${PROJECT_BINARY_DIR}")

OPTION(DEBUG_ALL "Include all debug code" OFF)

#
# Add subdirectories
#
add_subdirectory( src )
add_subdirectory( tests )

#
# configure a header file to pass cmake variables to source code
#
configure_file ( "${PROJECT_SOURCE_DIR}/LpmConfig.h.in" "${PROJECT_BINARY_DIR}/LpmConfig.h" @ONLY)

ENABLE_TESTING()
