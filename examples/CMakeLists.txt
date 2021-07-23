# Author: José D. Tascón-Vidarte <josedtascon@gmail.com>

include_directories(${CMAKE_INCLUDE_DIRECTORIES_BEFORE})
# include_directories(${PROJECT_SOURCE_DIR}/inc)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/cameraman20x16.tif
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/sinc_pad.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/scaled_rire_ct_001.nrrd
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/lenag0.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/lenag1.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/t0.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/t1.png
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/lpba_s01.nrrd
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/lpba_s02.nrrd
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

# ===========================================
#                  Build
# ===========================================
# Core
# Extra layer used during development
option(ENABLE_SAMPLES "Enable Samples" OFF) 
if(ENABLE_SAMPLES)

add_executable(image2d_arithmetic image2d_arithmetic.cpp)
target_link_libraries(image2d_arithmetic ${link_targets})

add_executable(image_read_write image_read_write.cpp)
target_link_libraries(image_read_write ${link_targets})

add_executable(image_statistics_input image_statistics_input.cpp)
target_link_libraries(image_statistics_input ${link_targets})

add_executable(image_fftconvolve_cpu image_fftconvolve_cpu.cpp)
target_link_libraries(image_fftconvolve_cpu ${link_targets})

add_executable(image_utils_cpu image_utils_cpu.cpp)
target_link_libraries(image_utils_cpu ${link_targets})

add_executable(vector_cpu_aritmethic vector_cpu_aritmethic.cpp)
target_link_libraries(vector_cpu_aritmethic ${link_targets})

# add_executable(vector_vcl_aritmethic vector_vcl_aritmethic.cpp)
# target_link_libraries(vector_vcl_aritmethic ${link_targets})

add_executable(image_gradient_kernel image_gradient_kernel.cpp)
target_link_libraries(image_gradient_kernel ${link_targets})

# # 2d test
# add_executable(image2d_basic image2d_basic.cpp)
# target_link_libraries(image2d_basic ${link_targets})

add_executable(image_assert image_assert.cpp)
target_link_libraries(image_assert ${link_targets})

add_executable(image_copy image_copy.cpp)
target_link_libraries(image_copy ${link_targets})

# add_executable(image2d_complex image2d_complex.cpp)
# target_link_libraries(image2d_complex ${link_targets})

add_executable(image_conditions_cpu image_conditions_cpu.cpp)
target_link_libraries(image_conditions_cpu ${link_targets})


add_executable(image_reduction_cpu image_reduction_cpu.cpp)
target_link_libraries(image_reduction_cpu ${link_targets})

add_executable(image_region image_region.cpp)
target_link_libraries(image_region ${link_targets})

add_executable(image_bounding_box image_bounding_box.cpp)
target_link_libraries(image_bounding_box ${link_targets})

# add_executable(image2d_fft image2d_fft.cpp)
# target_link_libraries(image2d_fft ${link_targets})

add_executable(image_fft_cpu image_fft_cpu.cpp)
target_link_libraries(image_fft_cpu ${link_targets})

add_executable(image_gradient_cpu image_gradient_cpu.cpp)
target_link_libraries(image_gradient_cpu ${link_targets})

add_executable(image2d_itk image2d_itk.cpp)
target_link_libraries(image2d_itk ${link_targets})

# add_executable(image2d_benchmarks image2d_benchmarks.cpp)
# target_link_libraries(image2d_benchmarks ${link_targets} "benchmark")

add_executable(image_meshgrid_cpu image_meshgrid_cpu.cpp)
target_link_libraries(image_meshgrid_cpu ${link_targets})

add_executable(affine_basic_cpu affine_basic_cpu.cpp)
target_link_libraries(affine_basic_cpu ${link_targets})

add_executable(dfield_basic_cpu dfield_basic_cpu.cpp)
target_link_libraries(dfield_basic_cpu ${link_targets})

add_executable(dfield_read_write_cpu dfield_read_write_cpu.cpp)
target_link_libraries(dfield_read_write_cpu ${link_targets})

add_executable(dfield_transform_input dfield_transform_input.cpp)
target_link_libraries(dfield_transform_input ${link_targets})

add_executable(jacobian_dfield_transform_input jacobian_dfield_transform_input.cpp)
target_link_libraries(jacobian_dfield_transform_input ${link_targets})

add_executable(transform_operations_cpu transform_operations_cpu.cpp)
target_link_libraries(transform_operations_cpu ${link_targets})

add_executable(transform_affine_read_write transform_affine_read_write.cpp)
target_link_libraries(transform_affine_read_write ${link_targets})

# add_executable(interpolator_basic interpolator_basic.cpp)
# target_link_libraries(interpolator_basic ${link_targets})

add_executable(image_interpolator_nearest_cpu image_interpolator_nearest_cpu.cpp)
target_link_libraries(image_interpolator_nearest_cpu ${link_targets})

add_executable(image_interpolator_linear_cpu image_interpolator_linear_cpu.cpp)
target_link_libraries(image_interpolator_linear_cpu ${link_targets})


# add_executable(image2d_interpolator_linear image2d_interpolator_linear.cpp)
# target_link_libraries(image2d_interpolator_linear ${link_targets})

# add_executable(transform_compare_itk transform_compare_itk.cpp)
# target_link_libraries(transform_compare_itk ${link_targets})

# add_executable(image2d_metric image2d_metric.cpp)
# target_link_libraries(image2d_metric ${link_targets})

add_executable(image_metric_cpu image_metric_cpu.cpp)
target_link_libraries(image_metric_cpu ${link_targets})

add_executable(image_resolution image_resolution.cpp)
target_link_libraries(image_resolution ${link_targets})

add_executable(image_resolution_input2d image_resolution_input2d.cpp)
target_link_libraries(image_resolution_input2d ${link_targets})

add_executable(image_resolution_input3d image_resolution_input3d.cpp)
target_link_libraries(image_resolution_input3d ${link_targets})


add_executable(image_optimizer_cpu image_optimizer_cpu.cpp)
target_link_libraries(image_optimizer_cpu ${link_targets})

# if(ENABLE_SAMPLES)

add_executable(image_optimizer_demons_cpu image_optimizer_demons_cpu.cpp)
target_link_libraries(image_optimizer_demons_cpu ${link_targets})

add_executable(image_regularizer image_regularizer.cpp)
target_link_libraries(image_regularizer ${link_targets})

add_executable(gaussian_filter gaussian_filter.cpp)
target_link_libraries(gaussian_filter ${link_targets})

add_executable(registration_demons_cpu registration_demons_cpu.cpp)
target_link_libraries(registration_demons_cpu ${link_targets})

add_executable(registration_demons_input registration_demons_input.cpp)
target_link_libraries(registration_demons_input ${link_targets})

add_executable(registration_demonsloi_input registration_demonsloi_input.cpp)
target_link_libraries(registration_demonsloi_input ${link_targets})

add_executable(registration_lddmm_prototype registration_lddmm_prototype.cpp)
target_link_libraries(registration_lddmm_prototype ${link_targets})

add_executable(registration_lddmm_input registration_lddmm_input.cpp)
target_link_libraries(registration_lddmm_input ${link_targets})

# endif()

if (IMART_WITH_OPENCL)

    # add_executable(gpu_image2d_basic gpu_image2d_basic.cpp)
    # target_link_libraries(gpu_image2d_basic ${link_targets} )

    # add_executable(gpu_image2d_arithmetic gpu_image2d_arithmetic.cpp)
    # target_link_libraries(gpu_image2d_arithmetic ${link_targets} )

    add_executable(opencl_object opencl_object.cpp)
    target_link_libraries(opencl_object ${link_targets} )

    add_executable(vector_opencl_aritmethic vector_opencl_aritmethic.cpp)
    target_link_libraries(vector_opencl_aritmethic ${link_targets} )

    add_executable(image2d_arithmetic_opencl image2d_arithmetic_opencl.cpp)
    target_link_libraries(image2d_arithmetic_opencl ${link_targets} )

    add_executable(image_conditions_opencl image_conditions_opencl.cpp)
    target_link_libraries(image_conditions_opencl ${link_targets} )

    add_executable(image_utils_opencl image_utils_opencl.cpp)
    target_link_libraries(image_utils_opencl ${link_targets} )

    add_executable(image_reduction_opencl image_reduction_opencl.cpp)
    target_link_libraries(image_reduction_opencl ${link_targets} )

    if (IMART_WITH_CLFFT)
    add_executable(image_fft_opencl image_fft_opencl.cpp)
    target_link_libraries(image_fft_opencl ${link_targets} )

    add_executable(image_gradient_opencl_input image_gradient_opencl_input.cpp)
    target_link_libraries(image_gradient_opencl_input ${link_targets} )
    endif()

    add_executable(image_gradient_kernel_opencl image_gradient_kernel_opencl.cpp)
    target_link_libraries(image_gradient_kernel_opencl ${link_targets} )

    add_executable(image_gradient_opencl image_gradient_opencl.cpp)
    target_link_libraries(image_gradient_opencl ${link_targets} )

    add_executable(image_meshgrid_opencl image_meshgrid_opencl.cpp)
    target_link_libraries(image_meshgrid_opencl ${link_targets} )

    add_executable(affine_basic_opencl affine_basic_opencl.cpp)
    target_link_libraries(affine_basic_opencl ${link_targets} )

    add_executable(dfield_basic_opencl dfield_basic_opencl.cpp)
    target_link_libraries(dfield_basic_opencl ${link_targets} )

    add_executable(image_interpolator_nearest_opencl image_interpolator_nearest_opencl.cpp)
    target_link_libraries(image_interpolator_nearest_opencl ${link_targets} )

    add_executable(image_interpolator_linear_opencl image_interpolator_linear_opencl.cpp)
    target_link_libraries(image_interpolator_linear_opencl ${link_targets} )

    add_executable(image_metric_opencl image_metric_opencl.cpp)
    target_link_libraries(image_metric_opencl ${link_targets} )

    add_executable(image_optimizer_opencl image_optimizer_opencl.cpp)
    target_link_libraries(image_optimizer_opencl ${link_targets} )

    add_executable(image_optimizer_demons_opencl image_optimizer_demons_opencl.cpp)
    target_link_libraries(image_optimizer_demons_opencl ${link_targets} )

    add_executable(registration_demons_opencl registration_demons_opencl.cpp)
    target_link_libraries(registration_demons_opencl ${link_targets} )

    add_executable(registration_opencl_demons_input registration_opencl_demons_input.cpp)
    target_link_libraries(registration_opencl_demons_input ${link_targets} )
endif()

if (IMART_WITH_CUDA)
    add_executable(vector_cuda_aritmethic vector_cuda_aritmethic.cpp)
    target_link_libraries(vector_cuda_aritmethic ${link_targets} )

    add_executable(image_utils_cuda image_utils_cuda.cpp)
    target_link_libraries(image_utils_cuda ${link_targets} )

    add_executable(image_reduction_cuda image_reduction_cuda.cpp)
    target_link_libraries(image_reduction_cuda ${link_targets} )

    add_executable(image_gradient_cuda image_gradient_cuda.cpp)
    target_link_libraries(image_gradient_cuda ${link_targets} )

    add_executable(image_fft_cuda image_fft_cuda.cpp)
    target_link_libraries(image_fft_cuda ${link_targets} )

    add_executable(image_meshgrid_cuda image_meshgrid_cuda.cpp)
    target_link_libraries(image_meshgrid_cuda ${link_targets} )

    add_executable(affine_basic_cuda affine_basic_cuda.cpp)
    target_link_libraries(affine_basic_cuda ${link_targets} )

    add_executable(dfield_basic_cuda dfield_basic_cuda.cpp)
    target_link_libraries(dfield_basic_cuda ${link_targets} )

    add_executable(image_interpolator_nearest_cuda image_interpolator_nearest_cuda.cpp)
    target_link_libraries(image_interpolator_nearest_cuda ${link_targets} )

    add_executable(image_interpolator_linear_cuda image_interpolator_linear_cuda.cpp)
    target_link_libraries(image_interpolator_linear_cuda ${link_targets} )

    add_executable(image_optimizer_cuda image_optimizer_cuda.cpp)
    target_link_libraries(image_optimizer_cuda ${link_targets} )

    add_executable(image_optimizer_demons_cuda image_optimizer_demons_cuda.cpp)
    target_link_libraries(image_optimizer_demons_cuda ${link_targets} )

    add_executable(registration_demons_cuda registration_demons_cuda.cpp)
    target_link_libraries(registration_demons_cuda ${link_targets} )

    add_executable(registration_cuda_demons_input registration_cuda_demons_input.cpp)
    target_link_libraries(registration_cuda_demons_input ${link_targets} )

    # add_executable(registration_cuda_lddmm_prototype registration_cuda_lddmm_prototype.cpp)
    # target_link_libraries(registration_cuda_lddmm_prototype ${link_targets} )
endif()

# # 3d test
# add_executable(image3d_basic image3d_basic.cpp)
# target_link_libraries(image3d_basic ${link_targets})

# add_executable(image3d_itk image3d_itk.cpp)
# target_link_libraries(image3d_itk ${link_targets})


if (IMART_WITH_OPENCL AND IMART_WITH_CUDA)
    # CHECK if OPENCL and CUDA
    add_executable(affine_ilinear_cpu_gpu_itk affine_ilinear_cpu_gpu_itk.cpp)
    target_link_libraries(affine_ilinear_cpu_gpu_itk ${link_targets}  )

    add_executable(scale_icubic_cpu_gpu_itk scale_icubic_cpu_gpu_itk.cpp)
    target_link_libraries(scale_icubic_cpu_gpu_itk ${link_targets}  )

    add_executable(transform_dfield_cpu_gpu transform_dfield_cpu_gpu.cpp)
    target_link_libraries(transform_dfield_cpu_gpu ${link_targets}  )
endif()

# Real time
add_executable(viewer viewer.cpp)
target_link_libraries(viewer ${link_targets})

add_executable(tracking_viewer tracking_viewer.cpp)
target_link_libraries(tracking_viewer ${link_targets})

endif()

# add_executable(registration_horn_schunck_input registration_horn_schunck_input.cpp)
# target_link_libraries(registration_horn_schunck_input ${link_targets})

add_executable(tracking_simulation_registration tracking_simulation_registration.cpp)
target_link_libraries(tracking_simulation_registration ${link_targets})

add_executable(tracking_viewer_multiple_organs tracking_viewer_multiple_organs.cpp)
target_link_libraries(tracking_viewer_multiple_organs ${link_targets})

add_executable(tracking_template_matching tracking_template_matching.cpp)
target_link_libraries(tracking_template_matching ${link_targets})

add_executable(tracking_template_matching_and_registration tracking_template_matching_and_registration.cpp)
target_link_libraries(tracking_template_matching_and_registration ${link_targets})

