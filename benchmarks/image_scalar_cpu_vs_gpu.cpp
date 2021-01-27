/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-04 10:35:41
*/

#include <iostream>
#include <benchmark/benchmark.h>

// local headers
#include "../src/image.h"

using namespace imart;

// Function to be timed
static void scalar_add_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cpu<type> img1(N,N);
    image_cpu<type> img2(N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 + scalar;
    };
};

static void scalar_add_image_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_gpu<type> img1(N,N);
    image_gpu<type> img2(N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 + scalar;
    };
};

static void scalar_add_image_cuda_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cuda<type> img1(N,N);
    image_cuda<type> img2(N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 + scalar;
    };
};

static void scalar_add_image_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cpu<type> img1(N,N,N);
    image_cpu<type> img2(N,N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 + scalar;
    };
};

static void scalar_add_image_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_gpu<type> img1(N,N,N);
    image_gpu<type> img2(N,N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 + scalar;
    };
};

static void scalar_add_image_cuda_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cuda<type> img1(N,N,N);
    image_cuda<type> img2(N,N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 + scalar;
    };
};


// Register the function as a benchmark
BENCHMARK(scalar_add_image_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(scalar_add_image_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(scalar_add_image_cuda_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(scalar_add_image_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(scalar_add_image_gpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(scalar_add_image_cuda_3d)->RangeMultiplier(10)->Range(10, 400);


// Run the benchmark
BENCHMARK_MAIN();
