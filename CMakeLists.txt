# Author: José D. Tascón-Vidarte <josedtascon@gmail.com>

cmake_minimum_required(VERSION 3.8)
project(imart C CXX)

set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)
set(VERSION_PATCH 3)
set(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Path where all your libraries are
set(CMAKE_LIBRARY_PATH "/usr/lib")

# ===========================================
#               Find Libraries
# ===========================================
# Find ITK
find_package( ITK REQUIRED )
include(${ITK_USE_FILE})
# Find OpenMP
# find_package( OpenMP REQUIRED )
# Find OpenCL.
find_package( OpenCL REQUIRED )


set(link_targets ${ITK_LIBRARIES} OpenCL::OpenCL fftw3  gtest)
# set(link_targets OpenMP::OpenMP_CXX )

# ===========================================
#               Compiler Flags
# ===========================================
set (CMAKE_CXX_STANDARD 17)



# ===========================================
#                  Options
# ===========================================
option(BUILD_TESTS "Build the tests" ON)

if(BUILD_TESTS)
   add_subdirectory(tests)
endif()


# ===========================================
#                  Build
# ===========================================
# Core includes
include_directories(src src/utils)
add_subdirectory(src)
# subdirs(src)

# add_executable( image_register src )
# target_link_libraries( image_register src )

# Link 
# ITK = ${ITK_LIBRARIES}  
# OpenCL = ${OPENCL_LIBRARIES} OpenCL::OpenCL
# OpenMP = OpenMP::OpenMP_CXX

# ===========================================
#                  Install
# ===========================================

# to include image utils in the future
# install(TARGETS image_info DESTINATION bin)
# install(TARGETS image_register DESTINATION bin)

# this ones require vtk
# install(TARGETS image_show2d DESTINATION bin)
# install(TARGETS image_show3d DESTINATION bin)