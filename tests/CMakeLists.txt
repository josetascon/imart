include_directories(${CMAKE_INCLUDE_DIRECTORIES_BEFORE})
# include_directories(${PROJECT_SOURCE_DIR}/inc)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/cameraman20x16.tif
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

# add_library(tests image2d_basic.cpp image2d_itk.cpp)

# 2d test
add_executable(image2d_basic image2d_basic.cpp)
target_link_libraries(image2d_basic src ${link_targets})

add_executable(image2d_itk image2d_itk.cpp)
target_link_libraries(image2d_itk src ${link_targets})

add_executable(image2d_arithmetic image2d_arithmetic.cpp)
target_link_libraries(image2d_arithmetic src ${link_targets})

add_executable(image2d_benchmarks image2d_benchmarks.cpp)
target_link_libraries(image2d_benchmarks src ${link_targets} "benchmark")

add_executable(image2d_meshgrid image2d_meshgrid.cpp)
target_link_libraries(image2d_meshgrid src ${link_targets})

add_executable(transform2d_basic transform2d_basic.cpp)
target_link_libraries(transform2d_basic src ${link_targets})

add_executable(image2d_interpolate image2d_interpolate.cpp)
target_link_libraries(image2d_interpolate src ${link_targets})

# # 3d test
add_executable(image3d_basic image3d_basic.cpp)
target_link_libraries(image3d_basic src ${link_targets})