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

add_executable(vector_vcl vector_vcl.cpp)
target_link_libraries(vector_vcl ${test_targets})

add_executable(vector_ocl vector_ocl.cpp)
target_link_libraries(vector_ocl ${test_targets})

add_executable(vector_cuda vector_cuda.cpp)
target_link_libraries(vector_cuda ${test_targets})


# Images
add_executable(image_cpu image_cpu.cpp)
target_link_libraries(image_cpu ${test_targets})



# Transformations
add_executable(dfield_cpu dfield_cpu.cpp)
target_link_libraries(dfield_cpu ${test_targets})

add_executable(dfield_cuda dfield_cuda.cpp)
target_link_libraries(dfield_cuda ${test_targets})