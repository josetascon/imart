/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-26 16:02:43
*/

#include <iostream>
#include <benchmark/benchmark.h>

// local headers
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"
#include "../src/ilinear.h"

using namespace imart;

// Function to be timed
static void bm_ilinear_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_opencl<type>::new_pointer(N,N);
    auto image1 = image_opencl<type>::new_pointer();
    image0->random();
    auto x0 = grid_opencl<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_opencl<type>::pointer params(new image_opencl<type>{1.1, 0.1, -0.2, 0.9, 1.3, 8.0});
    auto taffine = affine_opencl<type>::new_pointer(2, params);
    auto interp0 = ilinear_opencl<type>::new_pointer(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply(x0);
        image1 = interp0->apply(x1);
    };
    // x1.print_data();
};

static void bm_ilinear_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_opencl<type>::new_pointer(N,N,N);
    auto image1 = image_opencl<type>::new_pointer();
    image0->random();
    auto x0 = grid_opencl<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_opencl<type>::pointer params(new image_opencl<type>{1.1, 0.1, -0.2, 0.05, 1.2, 0.03, 0, -0.04, 1, 11.324, 201.4, 8.0});
    auto taffine = affine_opencl<type>::new_pointer(3, params);
    auto interp0 = ilinear_opencl<type>::new_pointer(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply(x0);
        image1 = interp0->apply(x1);
    };
    // x1.print_data();
};

// Register the function as a benchmark
// BENCHMARK(bm_ilinear_opencl_2d)->RangeMultiplier(10)->Range(10, 8000);
// BENCHMARK(bm_ilinear_cuda_2d)->RangeMultiplier(10)->Range(10, 8000);
// BENCHMARK(bm_ilinear_opencl_3d)->RangeMultiplier(10)->Range(10, 400);
// BENCHMARK(bm_ilinear_cuda_3d)->RangeMultiplier(10)->Range(10, 400);

static void CustomArguments2d(benchmark::internal::Benchmark* b)
{
    std::vector<int> list = {10,20,40,80,100,200,400,800,1000,2000,4000,8000};
    for(int i = 0; i < list.size(); i++)
        b->Args({list[i]});
};

static void CustomArguments3d(benchmark::internal::Benchmark* b)
{
    std::vector<int> list = {10,20,40,80,100,200,400,800};
    for(int i = 0; i < list.size(); i++)
        b->Args({list[i]});
};

BENCHMARK(bm_ilinear_opencl_2d)->Apply(CustomArguments2d);
BENCHMARK(bm_ilinear_opencl_3d)->Apply(CustomArguments3d);


// Run the benchmark
BENCHMARK_MAIN();