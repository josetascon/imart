# Author: José D. Tascón-Vidarte <josedtascon@gmail.com>
# ===========================================
#               CMAKE Project
# ===========================================
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(imart C CXX)
# project(imart LANGUAGES C CXX CUDA)

# set(VERSION_MAJOR 0)
# set(VERSION_MINOR 2)
# set(VERSION_PATCH 1)
# set(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Path where all your libraries are
set(CMAKE_LIBRARY_PATH "/usr/lib")

# ===========================================
#               Dependencies
# ===========================================
option(ENABLE_OPENMP "Enable OpenMP" ON)
option(ENABLE_OPENCL "Enable OpenCL" ON)
option(ENABLE_CUDA "Enable CUDA" ON)

option(ENABLE_FFTW "Enable fftw3" ON)
option(ENABLE_CLFFTW "Enable clfft" ON)


# ===========================================
#               Compiler Flags
# ===========================================
# c++ standard
set (CMAKE_CXX_STANDARD 17)
# set (CMAKE_CUDA_STANDARD 17)
set (CMAKE_CXX_FLAGS "-export-dynamic")

# build type
set (CMAKE_BUILD_TYPE "Release")
# set (CMAKE_BUILD_TYPE "Debug")
# set (CMAKE_BUILD_TYPE "RelWithDebInfo")

if(ENABLE_OPENMP)
    find_package( OpenMP REQUIRED )
endif()



# ===========================================
#               Find Libraries
# ===========================================
# Find ITK
find_package( ITK REQUIRED )
include(${ITK_USE_FILE})

find_package( Boost COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )


if(ENABLE_OPENCL)
    # Find OpenCL.
    find_package( OpenCL REQUIRED )
endif()

# ===========================================
#                  Includes
# ===========================================
# Core includes
include_directories(src src/utils)

# ===========================================
#              Libs: imart_src
# ===========================================
add_subdirectory(src)
# subdirs(src)

# ===========================================
#              Libs: imart_cuda
# ===========================================
if(ENABLE_CUDA)
    enable_language(CUDA)
    set (CMAKE_CUDA_ARCHITECTURES 52 72)
    add_subdirectory(src/cuda)
endif()

# ===========================================
#                   Links
# ===========================================
# set(link_targets imart_src ${ITK_LIBRARIES})
set(link_targets imart_src ${ITK_LIBRARIES} ${Boost_LIBRARIES} dl)

if (ENABLE_OPENMP)
    set(link_targets ${link_targets} OpenMP::OpenMP_CXX)
endif()
if (ENABLE_OPENCL)
    set(link_targets ${link_targets} OpenCL::OpenCL)
endif()
if (ENABLE_CUDA)
    set(link_targets ${link_targets} imart_cuda)
endif()

if (ENABLE_FFTW)
    set(link_targets ${link_targets} fftw3)
endif()
if (ENABLE_CLFFTW)
    set(link_targets ${link_targets} clFFT)
endif()

# ===========================================
#               Build Options
# ===========================================
option(BUILD_TESTS "Build the tests" OFF)

if(BUILD_TESTS)
   add_subdirectory(tests)
endif()

option(BUILD_EXAMPLES "Build the examples" ON)
if(BUILD_EXAMPLES)
   add_subdirectory(examples)
endif()

option(BUILD_BENCHMARKS "Build the benchmarks" OFF)
if(BUILD_BENCHMARKS)
   add_subdirectory(benchmarks)
endif()

# ===========================================
#               Command Build
# ===========================================
# add_executable( image_register src )
# target_link_libraries( image_register src )

# ===========================================
#                  Install
# ===========================================

# to include image utils in the future
# install(TARGETS image_info DESTINATION bin)
# install(TARGETS image_register DESTINATION bin)

# this ones require vtk
# install(TARGETS image_show2d DESTINATION bin)
# install(TARGETS image_show3d DESTINATION bin)