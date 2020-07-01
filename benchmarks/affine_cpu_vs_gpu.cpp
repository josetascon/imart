/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-29 10:17:32
*/

#include <iostream>
#include <benchmark/benchmark.h>

// local headers
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"

using namespace imart;

// Function to be timed
static void bm_affine_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_cpu<type>::new_pointer(N,N);
    auto x0 = grid_cpu<type>::new_pointer(*image0);
    auto x1 = x0->mimic();
    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.9, 1.3, 8.0});
    auto taffine = affine<type>::new_pointer(2, params);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply<vector_cpu<type>>(x0);
    };
    // x1->print_data();
};

static void bm_affine_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_cpu<type>::new_pointer(N,N,N);
    auto x0 = grid_cpu<type>::new_pointer(*image0);
    auto x1 = x0->mimic();
    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.05, 1.2, 0.03, 0, -0.04, 1, 11.324, 201.4, 8.0});
    auto taffine = affine<type>::new_pointer(3, params);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply<vector_cpu<type>>(x0);
    };
    // x1->print_data();
};

static void bm_affine_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_gpu<type>::new_pointer(N,N);
    auto x0 = grid_gpu<type>::new_pointer(*image0);
    auto x1 = x0->mimic();
    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.9, 1.3, 8.0});
    auto taffine = affine<type>::new_pointer(2, params);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply<vector_ocl<type>>(x0);
    };
    // x1->print_data();
};

static void bm_affine_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_gpu<type>::new_pointer(N,N,N);
    auto x0 = grid_gpu<type>::new_pointer(*image0);
    auto x1 = x0->mimic();
    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.05, 1.2, 0.03, 0, -0.04, 1, 11.324, 201.4, 8.0});
    auto taffine = affine<type>::new_pointer(3, params);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply<vector_ocl<type>>(x0);
    };
    // x1->print_data();
};

// Register the function as a benchmark
BENCHMARK(bm_affine_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_affine_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_affine_cpu_3d)->RangeMultiplier(10)->Range(10, 300);
BENCHMARK(bm_affine_gpu_3d)->RangeMultiplier(10)->Range(10, 300);



// Run the benchmark
BENCHMARK_MAIN();
