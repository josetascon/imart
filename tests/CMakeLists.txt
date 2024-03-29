# Author: José D. Tascón-Vidarte <josedtascon@gmail.com>

include_directories(${CMAKE_INCLUDE_DIRECTORIES_BEFORE})

set(test_targets ${link_targets} gtest)

# ===========================================
#                  Build
# ===========================================
# Core

# Vectors
add_executable(vector_cpu vector_cpu.cpp)
target_link_libraries(vector_cpu ${test_targets})

if (IMART_WITH_OPENCL)
    add_executable(vector_opencl vector_opencl.cpp)
    target_link_libraries(vector_opencl ${test_targets})
endif()

if (IMART_WITH_CUDA)
    add_executable(vector_cuda vector_cuda.cpp)
    target_link_libraries(vector_cuda ${test_targets})
endif()

# add_executable(vector_vcl vector_vcl.cpp)
# target_link_libraries(vector_vcl ${test_targets})

# Images
add_executable(image_cpu image_cpu.cpp)
target_link_libraries(image_cpu ${test_targets})


# Transformations
add_executable(dfield_cpu dfield_cpu.cpp)
target_link_libraries(dfield_cpu ${test_targets})

if (IMART_WITH_CUDA)
    add_executable(dfield_cuda dfield_cuda.cpp)
    target_link_libraries(dfield_cuda ${test_targets})
endif()