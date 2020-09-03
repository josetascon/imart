# Author: José D. Tascón-Vidarte <josedtascon@gmail.com>

include_directories(${CMAKE_INCLUDE_DIRECTORIES_BEFORE})
# include_directories(${PROJECT_SOURCE_DIR}/inc)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/cameraman20x16.tif
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/sinc_pad.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/scaled_rire_ct_001.nrrd
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/lenag1.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/lenag2.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

# ===========================================
#                  Build
# ===========================================
# Core
add_executable(opencl_object opencl_object.cpp)
target_link_libraries(opencl_object ${link_targets})

add_executable(vector_cpu_aritmethic vector_cpu_aritmethic.cpp)
target_link_libraries(vector_cpu_aritmethic ${link_targets})

add_executable(vector_vcl_aritmethic vector_vcl_aritmethic.cpp)
target_link_libraries(vector_vcl_aritmethic ${link_targets})

add_executable(vector_ocl_aritmethic vector_ocl_aritmethic.cpp)
target_link_libraries(vector_ocl_aritmethic ${link_targets})

add_executable(image_gradient_kernel image_gradient_kernel.cpp)
target_link_libraries(image_gradient_kernel ${link_targets})

# # 2d test
# add_executable(image2d_basic image2d_basic.cpp)
# target_link_libraries(image2d_basic ${link_targets})

# add_executable(image2d_arithmetic image2d_arithmetic.cpp)
# target_link_libraries(image2d_arithmetic ${link_targets})

# add_executable(image2d_arithmetic_gpu image2d_arithmetic_gpu.cpp)
# target_link_libraries(image2d_arithmetic_gpu ${link_targets})

# add_executable(image2d_complex image2d_complex.cpp)
# target_link_libraries(image2d_complex ${link_targets})

add_executable(image_utils_cpu image_utils_cpu.cpp)
target_link_libraries(image_utils_cpu ${link_targets})

add_executable(image_utils_gpu image_utils_gpu.cpp)
target_link_libraries(image_utils_gpu ${link_targets})

add_executable(image_reduction_cpu image_reduction_cpu.cpp)
target_link_libraries(image_reduction_cpu ${link_targets})

add_executable(image_reduction_gpu image_reduction_gpu.cpp)
target_link_libraries(image_reduction_gpu ${link_targets})

# add_executable(image2d_fft image2d_fft.cpp)
# target_link_libraries(image2d_fft ${link_targets})

add_executable(image_fft_cpu image_fft_cpu.cpp)
target_link_libraries(image_fft_cpu ${link_targets})

add_executable(image_fft_gpu image_fft_gpu.cpp)
target_link_libraries(image_fft_gpu ${link_targets})

add_executable(image_gradient_cpu image_gradient_cpu.cpp)
target_link_libraries(image_gradient_cpu ${link_targets})

add_executable(image_gradient_gpu image_gradient_gpu.cpp)
target_link_libraries(image_gradient_gpu ${link_targets})

add_executable(image_gradient_gpu_input image_gradient_gpu_input.cpp)
target_link_libraries(image_gradient_gpu_input ${link_targets})

add_executable(image2d_itk image2d_itk.cpp)
target_link_libraries(image2d_itk ${link_targets})

# add_executable(image2d_benchmarks image2d_benchmarks.cpp)
# target_link_libraries(image2d_benchmarks ${link_targets} "benchmark")

add_executable(image_meshgrid_cpu image_meshgrid_cpu.cpp)
target_link_libraries(image_meshgrid_cpu ${link_targets})

add_executable(image_meshgrid_gpu image_meshgrid_gpu.cpp)
target_link_libraries(image_meshgrid_gpu ${link_targets})

add_executable(affine_basic_cpu affine_basic_cpu.cpp)
target_link_libraries(affine_basic_cpu ${link_targets})

add_executable(affine_basic_gpu affine_basic_gpu.cpp)
target_link_libraries(affine_basic_gpu ${link_targets})

add_executable(dfield_basic_cpu dfield_basic_cpu.cpp)
target_link_libraries(dfield_basic_cpu ${link_targets})

add_executable(dfield_basic_gpu dfield_basic_gpu.cpp)
target_link_libraries(dfield_basic_gpu ${link_targets})

add_executable(transform_operations_cpu transform_operations_cpu.cpp)
target_link_libraries(transform_operations_cpu ${link_targets})

# add_executable(interpolator_basic interpolator_basic.cpp)
# target_link_libraries(interpolator_basic ${link_targets})

add_executable(image_interpolator_nearest_cpu image_interpolator_nearest_cpu.cpp)
target_link_libraries(image_interpolator_nearest_cpu ${link_targets})

add_executable(image_interpolator_linear_cpu image_interpolator_linear_cpu.cpp)
target_link_libraries(image_interpolator_linear_cpu ${link_targets})

add_executable(image_interpolator_nearest_gpu image_interpolator_nearest_gpu.cpp)
target_link_libraries(image_interpolator_nearest_gpu ${link_targets})

add_executable(image_interpolator_linear_gpu image_interpolator_linear_gpu.cpp)
target_link_libraries(image_interpolator_linear_gpu ${link_targets})

add_executable(transform_affine_cpu_gpu transform_affine_cpu_gpu.cpp)
target_link_libraries(transform_affine_cpu_gpu ${link_targets})

add_executable(transform_dfield_cpu_gpu transform_dfield_cpu_gpu.cpp)
target_link_libraries(transform_dfield_cpu_gpu ${link_targets})

# add_executable(image2d_interpolator_linear image2d_interpolator_linear.cpp)
# target_link_libraries(image2d_interpolator_linear ${link_targets})

# add_executable(image_copy image_copy.cpp)
# target_link_libraries(image_copy ${link_targets})

# add_executable(transform_compare_itk transform_compare_itk.cpp)
# target_link_libraries(transform_compare_itk ${link_targets})

# add_executable(image2d_metric image2d_metric.cpp)
# target_link_libraries(image2d_metric ${link_targets})

add_executable(image_metric_cpu image_metric_cpu.cpp)
target_link_libraries(image_metric_cpu ${link_targets})

add_executable(image_metric_gpu image_metric_gpu.cpp)
target_link_libraries(image_metric_gpu ${link_targets})

add_executable(image_resolution image_resolution.cpp)
target_link_libraries(image_resolution ${link_targets})

add_executable(image_resolution_input2d image_resolution_input2d.cpp)
target_link_libraries(image_resolution_input2d ${link_targets})

add_executable(image_resolution_input3d image_resolution_input3d.cpp)
target_link_libraries(image_resolution_input3d ${link_targets})

add_executable(image_optimizer_cpu image_optimizer_cpu.cpp)
target_link_libraries(image_optimizer_cpu ${link_targets})

add_executable(image_optimizer_gpu image_optimizer_gpu.cpp)
target_link_libraries(image_optimizer_gpu ${link_targets})

add_executable(image_optimizer_demons_cpu image_optimizer_demons_cpu.cpp)
target_link_libraries(image_optimizer_demons_cpu ${link_targets})

add_executable(image_optimizer_demons_gpu image_optimizer_demons_gpu.cpp)
target_link_libraries(image_optimizer_demons_gpu ${link_targets})

add_executable(gaussian_filter gaussian_filter.cpp)
target_link_libraries(gaussian_filter ${link_targets})

add_executable(registration_demons_cpu registration_demons_cpu.cpp)
target_link_libraries(registration_demons_cpu ${link_targets})

add_executable(registration_demons_gpu registration_demons_gpu.cpp)
target_link_libraries(registration_demons_gpu ${link_targets})


if (BUILD_CUDA)
    add_executable(vector_cuda_aritmethic vector_cuda_aritmethic.cpp)
    target_link_libraries(vector_cuda_aritmethic ${link_targets})

    add_executable(image_utils_cuda image_utils_cuda.cpp)
    target_link_libraries(image_utils_cuda ${link_targets})

    add_executable(image_reduction_cuda image_reduction_cuda.cpp)
    target_link_libraries(image_reduction_cuda ${link_targets})

    add_executable(image_gradient_cuda image_gradient_cuda.cpp)
    target_link_libraries(image_gradient_cuda ${link_targets})

    add_executable(image_meshgrid_cuda image_meshgrid_cuda.cpp)
    target_link_libraries(image_meshgrid_cuda ${link_targets})

    add_executable(affine_basic_cuda affine_basic_cuda.cpp)
    target_link_libraries(affine_basic_cuda ${link_targets})

    add_executable(dfield_basic_cuda dfield_basic_cuda.cpp)
    target_link_libraries(dfield_basic_cuda ${link_targets})

    add_executable(image_interpolator_nearest_cuda image_interpolator_nearest_cuda.cpp)
    target_link_libraries(image_interpolator_nearest_cuda ${link_targets})

    add_executable(image_interpolator_linear_cuda image_interpolator_linear_cuda.cpp)
    target_link_libraries(image_interpolator_linear_cuda ${link_targets})

    add_executable(image_optimizer_cuda image_optimizer_cuda.cpp)
    target_link_libraries(image_optimizer_cuda ${link_targets})

    add_executable(image_optimizer_demons_cuda image_optimizer_demons_cuda.cpp)
    target_link_libraries(image_optimizer_demons_cuda ${link_targets})

    add_executable(registration_demons_cuda registration_demons_cuda.cpp)
target_link_libraries(registration_demons_cuda ${link_targets})
endif()

# # 3d test
# add_executable(image3d_basic image3d_basic.cpp)
# target_link_libraries(image3d_basic ${link_targets})

# add_executable(image3d_itk image3d_itk.cpp)
# target_link_libraries(image3d_itk ${link_targets})


# # # GPU test
# add_executable(gpu_image2d_basic gpu_image2d_basic.cpp)
# target_link_libraries(gpu_image2d_basic ${link_targets})

# add_executable(gpu_image2d_arithmetic gpu_image2d_arithmetic.cpp)
# target_link_libraries(gpu_image2d_arithmetic ${link_targets})