
add_library(imart_cuda interface.cu interface.cuh kernels.cu kernels.cuh)

set_target_properties(imart_cuda PROPERTIES 
                              CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET imart_cuda PROPERTY CUDA_ARCHITECTURES 52 72)
# set_property(TARGET imart_cuda PROPERTY CUDA_STANDARD 17)

target_link_libraries(imart_cuda PUBLIC)