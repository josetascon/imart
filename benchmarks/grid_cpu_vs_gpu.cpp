/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-18 06:56:08
*/

#include <iostream>
#include <benchmark/benchmark.h>

// local headers
#include "../src/grid.h"

using namespace imart;

// Function to be timed
static void bm_grid_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cpu<type> image0(N,N);
    grid_cpu<type> x0(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x0.meshgrid();
    };
    // x0.print_data();
};

static void bm_grid_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cpu<type> image0(N,N,N);
    grid_cpu<type> x0(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x0.meshgrid();
    };
    // x0.print_data();
};

static void bm_grid_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);

    image_gpu<type> image0(N,N);
    grid_gpu<type> x0(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x0.meshgrid();
    };
    // x0.print_data();
};

static void bm_grid_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);

    image_gpu<type> image0(N,N,N);
    grid_gpu<type> x0(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x0.meshgrid();
    };
    // x0.print_data();
};

// Register the function as a benchmark
BENCHMARK(bm_grid_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_grid_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_grid_cpu_3d)->RangeMultiplier(10)->Range(10, 300);
BENCHMARK(bm_grid_gpu_3d)->RangeMultiplier(10)->Range(10, 300);



// Run the benchmark
BENCHMARK_MAIN();
