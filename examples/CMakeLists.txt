# Author: José D. Tascón-Vidarte <josedtascon@gmail.com>

include_directories(${CMAKE_INCLUDE_DIRECTORIES_BEFORE})
# include_directories(${PROJECT_SOURCE_DIR}/inc)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/cameraman20x16.tif
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/sinc_pad.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/scaled_rire_ct_001.nrrd
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

# ===========================================
#                  Build
# ===========================================
# Core
add_executable(opencl_object opencl_object.cpp)
target_link_libraries(opencl_object src ${link_targets})

add_executable(vector_cpu_aritmethic vector_cpu_aritmethic.cpp)
target_link_libraries(vector_cpu_aritmethic src ${link_targets})

add_executable(vector_vcl_aritmethic vector_vcl_aritmethic.cpp)
target_link_libraries(vector_vcl_aritmethic src ${link_targets})

add_executable(vector_ocl_aritmethic vector_ocl_aritmethic.cpp)
target_link_libraries(vector_ocl_aritmethic src ${link_targets})

add_executable(image_gradient_kernel image_gradient_kernel.cpp)
target_link_libraries(image_gradient_kernel src ${link_targets})

# # 2d test
# add_executable(image2d_basic image2d_basic.cpp)
# target_link_libraries(image2d_basic src ${link_targets})

# add_executable(image2d_arithmetic image2d_arithmetic.cpp)
# target_link_libraries(image2d_arithmetic src ${link_targets})

# add_executable(image2d_arithmetic_gpu image2d_arithmetic_gpu.cpp)
# target_link_libraries(image2d_arithmetic_gpu src ${link_targets})

# add_executable(image2d_complex image2d_complex.cpp)
# target_link_libraries(image2d_complex src ${link_targets})

add_executable(image_utils_cpu image_utils_cpu.cpp)
target_link_libraries(image_utils_cpu src ${link_targets})

add_executable(image_utils_gpu image_utils_gpu.cpp)
target_link_libraries(image_utils_gpu src ${link_targets})

add_executable(image_reduction_cpu image_reduction_cpu.cpp)
target_link_libraries(image_reduction_cpu src ${link_targets})

add_executable(image_reduction_gpu image_reduction_gpu.cpp)
target_link_libraries(image_reduction_gpu src ${link_targets})

# add_executable(image2d_fft image2d_fft.cpp)
# target_link_libraries(image2d_fft src ${link_targets})

add_executable(image_fft_cpu image_fft_cpu.cpp)
target_link_libraries(image_fft_cpu src ${link_targets})

add_executable(image_fft_gpu image_fft_gpu.cpp)
target_link_libraries(image_fft_gpu src ${link_targets})

add_executable(image_gradient_cpu image_gradient_cpu.cpp)
target_link_libraries(image_gradient_cpu src ${link_targets})

add_executable(image_gradient_gpu image_gradient_gpu.cpp)
target_link_libraries(image_gradient_gpu src ${link_targets})

add_executable(image_gradient_gpu_input image_gradient_gpu_input.cpp)
target_link_libraries(image_gradient_gpu_input src ${link_targets})

add_executable(image2d_itk image2d_itk.cpp)
target_link_libraries(image2d_itk src ${link_targets})

# add_executable(image2d_benchmarks image2d_benchmarks.cpp)
# target_link_libraries(image2d_benchmarks src ${link_targets} "benchmark")

add_executable(image_meshgrid_cpu image_meshgrid_cpu.cpp)
target_link_libraries(image_meshgrid_cpu src ${link_targets})

add_executable(image_meshgrid_gpu image_meshgrid_gpu.cpp)
target_link_libraries(image_meshgrid_gpu src ${link_targets})

# add_executable(transform2d_basic transform2d_basic.cpp)
# target_link_libraries(transform2d_basic src ${link_targets})

# add_executable(transform2d_basic_gpu transform2d_basic_gpu.cpp)
# target_link_libraries(transform2d_basic_gpu src ${link_targets})

# add_executable(interpolator_basic interpolator_basic.cpp)
# target_link_libraries(interpolator_basic src ${link_targets})

add_executable(image_interpolator_nearest_cpu image_interpolator_nearest_cpu.cpp)
target_link_libraries(image_interpolator_nearest_cpu src ${link_targets})

add_executable(image_interpolator_linear_cpu image_interpolator_linear_cpu.cpp)
target_link_libraries(image_interpolator_linear_cpu src ${link_targets})

add_executable(image_interpolator_nearest_gpu image_interpolator_nearest_gpu.cpp)
target_link_libraries(image_interpolator_nearest_gpu src ${link_targets})

add_executable(image_interpolator_linear_gpu image_interpolator_linear_gpu.cpp)
target_link_libraries(image_interpolator_linear_gpu src ${link_targets})

add_executable(transform_image_cpu_gpu transform_image_cpu_gpu.cpp)
target_link_libraries(transform_image_cpu_gpu src ${link_targets})

# add_executable(image2d_interpolator_linear image2d_interpolator_linear.cpp)
# target_link_libraries(image2d_interpolator_linear src ${link_targets})

# add_executable(image_copy image_copy.cpp)
# target_link_libraries(image_copy src ${link_targets})

# add_executable(transform_compare_itk transform_compare_itk.cpp)
# target_link_libraries(transform_compare_itk src ${link_targets})

# add_executable(image2d_metric image2d_metric.cpp)
# target_link_libraries(image2d_metric src ${link_targets})

add_executable(image_metric_cpu image_metric_cpu.cpp)
target_link_libraries(image_metric_cpu src ${link_targets})

add_executable(image_metric_gpu image_metric_gpu.cpp)
target_link_libraries(image_metric_gpu src ${link_targets})

add_executable(image_resolution image_resolution.cpp)
target_link_libraries(image_resolution src ${link_targets})

add_executable(image_resolution_input2d image_resolution_input2d.cpp)
target_link_libraries(image_resolution_input2d src ${link_targets})

add_executable(image_resolution_input3d image_resolution_input3d.cpp)
target_link_libraries(image_resolution_input3d src ${link_targets})

add_executable(image_optimizer_cpu image_optimizer_cpu.cpp)
target_link_libraries(image_optimizer_cpu src ${link_targets})

add_executable(image_optimizer_gpu image_optimizer_gpu.cpp)
target_link_libraries(image_optimizer_gpu src ${link_targets})

# # 3d test
# add_executable(image3d_basic image3d_basic.cpp)
# target_link_libraries(image3d_basic src ${link_targets})

# add_executable(image3d_itk image3d_itk.cpp)
# target_link_libraries(image3d_itk src ${link_targets})


# # # GPU test
# add_executable(gpu_image2d_basic gpu_image2d_basic.cpp)
# target_link_libraries(gpu_image2d_basic src ${link_targets})

# add_executable(gpu_image2d_arithmetic gpu_image2d_arithmetic.cpp)
# target_link_libraries(gpu_image2d_arithmetic src ${link_targets})