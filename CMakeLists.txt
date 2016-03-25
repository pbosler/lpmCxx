cmake_minimum_required ( VERSION 2.8 )

project ( LPM )
set (LPM_VERSION_MAJOR 2)
set (LPM_VERSION_MINOR 0)

set (CMAKE_BUILD_TYPE DEBUG)

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

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

#
# configure a header file to pass cmake variables to source code
#
configure_file ( "${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/lpmConfig.h")

include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SRC_DIR}/src" "${MPI_INCLUDE_PATH}" "${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src" "${MPI_INCLUDE_PATH}" "${PROJECT_BINARY_DIR}")

#
# Add subdirectories
#
add_subdirectory( src )
add_subdirectory( tests )


ENABLE_TESTING()
