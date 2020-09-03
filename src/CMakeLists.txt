add_library(imart_src main.cpp)

# set_target_properties(imart_src PROPERTIES 
#                               CUDA_SEPARABLE_COMPILATION ON)
# set_property(TARGET imart_src PROPERTY CUDA_ARCHITECTURES 52 72)

target_link_libraries(imart_src ${link_targets})