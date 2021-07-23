
add_library(imart_cuda interface.cu interface.cuh kernels.cu kernels.cuh)

# set_target_properties(imart_cuda PROPERTIES 
#                                 CUDA_SEPARABLE_COMPILATION ON
#                                 CUDA_ARCHITECTURES 86
#                                 CUDA_STANDARD 17)

set_target_properties(imart_cuda PROPERTIES 
                              CUDA_SEPARABLE_COMPILATION ON)
# set_property(TARGET imart_cuda PROPERTY CUDA_ARCHITECTURES 86)
# set_property(TARGET imart_cuda PROPERTY CUDA_STANDARD 17)

set_property(TARGET imart_cuda PROPERTY CUDA_ARCHITECTURES 52 72)

target_link_libraries(imart_cuda PUBLIC)

# target_link_libraries(imart_cuda PUBLIC CUDA::cudart)

# message([STATUS]"${link_targets}")

# target_link_libraries( imart_cuda PUBLIC ${link_targets} CUDA::cudart)

# target_include_directories( imart_cuda PUBLIC
#     ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# target_include_directories( imart_cuda PUBLIC
#                            ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )