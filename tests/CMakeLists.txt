include_directories(${CMAKE_INCLUDE_DIRECTORIES_BEFORE})
# include_directories(${PROJECT_SOURCE_DIR}/inc)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/cameraman20x16.tif
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images/)

add_library(tests image2d_basic.cpp image2d_itk.cpp)

add_executable(image2d_basic image2d_basic.cpp)
target_link_libraries(image2d_basic src tests ${link_targets})

add_executable(image2d_itk image2d_itk.cpp)
target_link_libraries(image2d_itk src tests ${link_targets})

add_executable(image2d_arithmetic image2d_arithmetic.cpp)
target_link_libraries(image2d_arithmetic src tests ${link_targets})

add_executable(image2d_benchmarks image2d_benchmarks.cpp)
target_link_libraries(image2d_benchmarks src tests ${link_targets} "benchmark")