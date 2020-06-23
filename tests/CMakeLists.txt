# Author: José D. Tascón-Vidarte <josedtascon@gmail.com>

include_directories(${CMAKE_INCLUDE_DIRECTORIES_BEFORE})

set(test_targets ${link_targets} gtest)

# ===========================================
#                  Build
# ===========================================
# Core
add_executable(vector_cpu vector_cpu.cpp)
target_link_libraries(vector_cpu src ${test_targets})

add_executable(vector_vcl vector_vcl.cpp)
target_link_libraries(vector_vcl src ${test_targets})

add_executable(vector_ocl vector_ocl.cpp)
target_link_libraries(vector_ocl src ${test_targets})

add_executable(image_cpu image_cpu.cpp)
target_link_libraries(image_cpu src ${test_targets})


