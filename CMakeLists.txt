cmake_minimum_required ( VERSION 3.1 )

project ( LPM )
set (LPM_VERSION_MAJOR 2)
set (LPM_VERSION_MINOR 0)

function (list2str list str)
  string (REPLACE ";" " " tmp "${list}")
  set (${str} ${tmp} PARENT_SCOPE)
endfunction ()

function (prc var)
  message ("${var}: ${${var}}")
endfunction ()

#
#   CMake 
#
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
set( CMAKE_CXX_FLAGS_DEBUG, "-g")
set( CMAKE_CXX_FLAGS_RELEASE, "-O3 -g -Wl,-stack_size -Wl,1000000")
option(USE_APPLE_CLANG "AppleClang compiler found. Disabling OpenMP." ON)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
set( CMAKE_CXX_FLAGS_DEBUG, "-g -fopenmp")
set( CMAKE_CXX_FLAGS_RELEASE, "-O3 -g -fopenmp")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
set( CMAKE_CXX_FLAGS_DEBUG, "-g -openmp")
set( CMAKE_CXX_FLAGS_RELEASE, "-O3 -g -openmp")
endif() 


SET (CMAKE_CXX_STANDARD 11) 

#
#   MPI
#
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
#
#   VTK
#
find_package(VTK REQUIRED HINTS $ENV{VTK_ROOT})
if (VTK_FOUND)
	OPTION(HAVE_VTK "Located vtk libraries" ON)
	Message("VTK DETAILS")
	Message("vtk_include_dirs = " ${VTK_INCLUDE_DIRS})
	Message("VTK_LIBRARIES = " ${VTK_LIBRARIES})
endif()
#
#   Kokkos
#
FIND_LIBRARY(Kokkos libkokkos.a REQUIRED HINTS ${Kokkos_ROOT}/lib $ENV{KOKKOS_ROOT}/lib)
if (Kokkos-NOTFOUND)
    Message(WARNING "WARNING: Kokkos library not found.")
endif()
FIND_PATH(KokkosInclude Kokkos_Core.hpp HINTS ${Kokkos_ROOT}/include $ENV{KOKKOS_ROOT}/include)
if (KokkosInclude-NOTFOUND)
    Message(WARNING "WARNING: Kokkos include directory not found.")
endif ()
if (Kokkos-NOTFOUND OR KokkosInclude-NOTFOUND) 
    Message(FATAL_ERROR "CMake could not find Kokkos.")
else ()
	OPTION(HAVE_KOKKOS "Located Kokkos library" ON)
	Message("KOKKOS DETAILS")
	Message("kokkos = " ${Kokkos})
	Message("kokkos-include = " ${KokkosInclude})
	include(${Kokkos_ROOT}/kokkos_generated_settings.cmake)
	list2str("${KOKKOS_CXXFLAGS}" Kokkos_CXX_FLAGS)
	list2str("${KOKKOS_LD_FLAGS}" Kokkos_LD_FLAGS)
	string (FIND "${KOKKOS_GMAKE_DEVICES}" "Cuda" cuda_str_pos)
	set (CUDA_BUILD FALSE)
	if (${cuda_str_pos} GREATER -1)
	    set (CUDA_BUILD TRUE)
	    execute_process(COMMAND ${CMAKE_CXX_COMPILER} "--nvcc-wrapper-show" 
	        RESULT_VARIABLE WRAPS_NVCC OUTPUT_VARIABLE WRAPS_NVCC_OUT)
        string (FIND ${WRAPS_NVCC_OUT} "nvcc" pos)
	endif()
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
