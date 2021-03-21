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
add_executable( memory_cpu2gpu memory_cpu2gpu.cpp )
target_link_libraries( memory_cpu2gpu ${bench_targets} )

add_executable( memory_gpu2cpu memory_gpu2cpu.cpp )
target_link_libraries( memory_gpu2cpu ${bench_targets} )

add_executable( image_copy_functions image_copy_functions.cpp )
target_link_libraries( image_copy_functions ${bench_targets} )

add_executable( image_scalar_functions image_scalar_functions.cpp )
target_link_libraries( image_scalar_functions ${bench_targets} )

add_executable( image_vector_operations image_vector_operations.cpp )
target_link_libraries( image_vector_operations ${bench_targets} )

add_executable( image_reduction_functions image_reduction_functions.cpp )
target_link_libraries( image_reduction_functions ${bench_targets} )

add_executable( image_fft image_fft.cpp )
target_link_libraries( image_fft ${bench_targets} )

add_executable( image_convolve image_convolve.cpp )
target_link_libraries( image_convolve ${bench_targets} )

add_executable( image_gradient_cpu_vs_gpu image_gradient_cpu_vs_gpu.cpp )
target_link_libraries( image_gradient_cpu_vs_gpu ${bench_targets} )

add_executable( grid_cpu_vs_gpu grid_cpu_vs_gpu.cpp )
target_link_libraries( grid_cpu_vs_gpu ${bench_targets} )

add_executable( affine_cpu_vs_gpu affine_cpu_vs_gpu.cpp )
target_link_libraries( affine_cpu_vs_gpu ${bench_targets} )

add_executable( interpolator_nearest_cpu_vs_gpu interpolator_nearest_cpu_vs_gpu.cpp )
target_link_libraries( interpolator_nearest_cpu_vs_gpu ${bench_targets} )

add_executable( interpolator_nearest_itk interpolator_nearest_itk.cpp )
target_link_libraries( interpolator_nearest_itk ${bench_targets} )

add_executable( interpolator_linear_cpu interpolator_linear_cpu.cpp )
target_link_libraries( interpolator_linear_cpu ${bench_targets} )

add_executable( interpolator_linear_itk interpolator_linear_itk.cpp )
target_link_libraries( interpolator_linear_itk ${bench_targets} )

add_executable( interpolator_linear_opencl interpolator_linear_opencl.cpp )
target_link_libraries( interpolator_linear_opencl ${bench_targets} )

add_executable( interpolator_linear_cuda interpolator_linear_cuda.cpp )
target_link_libraries( interpolator_linear_cuda ${bench_targets} )

add_executable( metric_ssd_cpu_vs_gpu metric_ssd_cpu_vs_gpu.cpp )
target_link_libraries( metric_ssd_cpu_vs_gpu ${bench_targets} )

add_executable( interpolator_nearest_cpu_gpu_cuda interpolator_nearest_cpu_gpu_cuda.cpp )
target_link_libraries( interpolator_nearest_cpu_gpu_cuda ${bench_targets} )

