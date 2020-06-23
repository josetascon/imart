# Author: José D. Tascón-Vidarte <josedtascon@gmail.com>

include_directories(${CMAKE_INCLUDE_DIRECTORIES_BEFORE})

# # ===========================================
# #               Find Libraries
# # ===========================================
# # Find benchmark
find_package( benchmark REQUIRED )

set(bench_targets ${link_targets} benchmark::benchmark)

# ===========================================
#                  Build
# ===========================================
add_executable( memory_cpu2gpu_vcl memory_cpu2gpu_vcl.cpp )
target_link_libraries( memory_cpu2gpu_vcl ${bench_targets} )

add_executable( memory_gpu2cpu_vcl memory_gpu2cpu_vcl.cpp )
target_link_libraries( memory_gpu2cpu_vcl ${bench_targets} )

add_executable( viennacl_vs_vector_vcl viennacl_vs_vector_vcl.cpp )
target_link_libraries( viennacl_vs_vector_vcl ${bench_targets} )

add_executable( vector_vcl_scalar vector_vcl_scalar.cpp )
target_link_libraries( vector_vcl_scalar ${bench_targets} )

add_executable( grid_cpu_vs_gpu grid_cpu_vs_gpu.cpp )
target_link_libraries( grid_cpu_vs_gpu ${bench_targets} )