cmake_minimum_required ( VERSION 3.1 )

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
SET (CMAKE_CXX_STANDARD 11) 

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

find_package(VTK REQUIRED HINTS $ENV{VTK_ROOT})
if (VTK_FOUND)
	OPTION(HAVE_VTK "Located vtk libraries" ON)
	Message("VTK DETAILS")
	Message("vtk_include_dirs = " ${VTK_INCLUDE_DIRS})
	Message("VTK_LIBRARIES = " ${VTK_LIBRARIES})
endif()

FIND_LIBRARY(Kokkos libkokkos.a HINTS $ENV{KOKKOS_ROOT}/lib $ENV{HOME}/kokkos-serial/install/lib)
if (Kokkos-NOTFOUND)
	Message(WARNING "WARNING: Kokkos library not found.")
else ()
	set(kokkos_lib ON)
endif()
FIND_PATH(KokkosInclude Kokkos_Core.hpp HINTS $ENV{KOKKOS_ROOT}/include)
if (KokkosInclude-NOTFOUND)
	Message(WARNING "WARNING: Kokkos include dir not found.")
else()
	set(kokkos_inc OFF)
endif()
if (Kokkos-NOTFOUND OR KokkosInclude-NOTFOUND) 
	Message(FATAL_ERROR "Kokkos not found." )
else()
	OPTION(HAVE_KOKKOS "Located Kokkos library" ON)
	Message("KOKKOS DETAILS")
	Message("kokkos = " ${Kokkos})
	Message("kokkos-include = " ${KokkosInclude})	
endif()

#
# Add subdirectories
#
add_subdirectory( src )
add_subdirectory( tests )
add_subdirectory( examples)

message("project binary dir = " ${PROJECT_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}" "${KokkosInclude}")


OPTION(DEBUG_ALL "Include all debug code" OFF)



#
# configure a header file to pass cmake variables to source code
#
configure_file ( "${CMAKE_SOURCE_DIR}/LpmConfig.h.in" "${PROJECT_BINARY_DIR}/LpmConfig.h" @ONLY)


ENABLE_TESTING()
