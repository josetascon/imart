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


# ===========================================
#               Options
# ===========================================
option(IMART_WITH_OPENMP "Enable OpenMP" ON)
option(IMART_WITH_OPENCL "Enable OpenCL" ON)
option(IMART_WITH_CUDA   "Enable CUDA"   ON)

option(IMART_WITH_FFTW   "Enable fftw3"  ON)
option(IMART_WITH_CLFFT  "Enable clfft"  ON)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(src/imart_config.h.in imart_config.h)


# ===========================================
#               Find Libraries
# ===========================================
# Find ITK
find_package( ITK REQUIRED )
include(${ITK_USE_FILE})

# Find VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package( Boost COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )


if(IMART_WITH_OPENMP)
    find_package( OpenMP REQUIRED)
endif()

if(IMART_WITH_OPENCL)
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
list(APPEND IMART_MODULES imart_src)
# subdirs(src)

# ===========================================
#              Libs: imart_cuda
# ===========================================
if(IMART_WITH_CUDA)
    enable_language(CUDA)
    set (CMAKE_CUDA_ARCHITECTURES 86) #52 72 )
    add_subdirectory(src/cuda)
endif()

# ===========================================
#                   Links
# ===========================================
# set(link_targets imart_src ${ITK_LIBRARIES})
# set(link_targets imart_src ${ITK_LIBRARIES} ${Boost_LIBRARIES} dl)
set(link_targets ${IMART_MODULES} ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${Boost_LIBRARIES} dl)

if (IMART_WITH_OPENMP)
    set(link_targets ${link_targets} OpenMP::OpenMP_CXX)
endif()

if (IMART_WITH_FFTW)
    set(link_targets ${link_targets} fftw3)
    if (IMART_WITH_OPENMP)
       set(link_targets ${link_targets} fftw3_omp)
    endif()
endif()

if (IMART_WITH_OPENCL)
    # list(APPEND opencl_targets OpenCL::OpenCL)
    set(link_targets ${link_targets} OpenCL::OpenCL)
endif()

if (IMART_WITH_CLFFT)
    # list(APPEND opencl_targets clFFT)
    set(link_targets ${link_targets} clFFT)
endif()

if (IMART_WITH_CUDA)
    # list(APPEND cuda_targets imart_cuda cufft)
    set(link_targets ${link_targets} imart_cuda cufft)
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
# Binaries
# install(TARGETS tracking_video_registration DESTINATION bin)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/imart_config.h" DESTINATION "/usr/include/imart")

# Headers
install(DIRECTORY "src/" # source directory
        DESTINATION "/usr/include/imart" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.h" # select header files
)