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
#cmake_policy(CMP0074 new)
prc(_CMAKE_TOOLCHAIN_PREFIX)
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
if (Kokkos_ROOT)
	include(${Kokkos_ROOT}/kokkos_generated_settings.cmake)
	set(KokkosInclude ${Kokkos_ROOT}/include)
	set(Kokkos ${Kokkos_ROOT}/lib/libkokkos.a)
	OPTION(HAVE_KOKKOS "Located Kokkos library" ON)
	Message("KOKKOS DETAILS")
	Message("kokkos = " ${Kokkos})
	Message("kokkos-include = " ${KokkosInclude})
	
	list2str("${KOKKOS_CXXFLAGS}" Kokkos_CXX_FLAGS)
	list2str("${KOKKOS_LD_FLAGS}" Kokkos_LD_FLAGS)
	prc(Kokkos_CXX_FLAGS)
	string (FIND "${KOKKOS_GMAKE_DEVICES}" "Cuda" cuda_str_pos)
	set (CUDA_BUILD FALSE)
	if (${cuda_str_pos} GREATER -1)
	    set (CUDA_BUILD TRUE)
	    execute_process(COMMAND ${CMAKE_CXX_COMPILER} "--nvcc-wrapper-show" 
	        RESULT_VARIABLE WRAPS_NVCC OUTPUT_VARIABLE WRAPS_NVCC_OUT)
        string (FIND ${WRAPS_NVCC_OUT} "nvcc" pos)
        prc(WRAPS_NVCC_OUT)
	endif()
else ()
	Message(FATAL_ERROR "Kokkos is required.  Set CMake variable Kokkos_ROOT.")
endif()
set(LPM_COMPILE_FLAGS "-g ${Kokkos_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
if (CUDA_BUILD)
	set(LPM_COMPILE_FLAGS "${LPM_COMPILE_FLAGS} --expt-relaxed-constexpr")
	get_property(cxx_features GLOBAL PROPERTY CMAKE_CXX_COMPILE_FEATURES)
	prc(CMAKE_CXX_COMPILE_FEATURES)
	foreach(feature ${cxx_features})
    message("${feature}")
endforeach()
	
endif()
set(LPM_LINK_FLAGS "-L${Kokkos_ROOT}/lib ${Kokkos_LD_FLAGS}")
prc(LPM_COMPILE_FLAGS)
prc(LPM_LINK_FLAGS)


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
