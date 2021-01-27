/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-04 19:02:53
*/

#include <iostream>
#include <benchmark/benchmark.h>

// local headers
#include "../src/image.h"

using namespace imart;

// Function to be timed
static void vector_add_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cpu<type> img1(N,N);
    image_cpu<type> img2(N,N);
    image_cpu<type> img3(N,N);
    img1.ones();
    img2.ones();
    img3.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img3 = img1 + img2;
    };
};

static void vector_add_image_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_gpu<type> img1(N,N);
    image_gpu<type> img2(N,N);
    image_gpu<type> img3(N,N);
    img1.ones();
    img2.ones();
    img3.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img3 = img1 + img2;
    };
};

static void vector_add_image_cuda_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cuda<type> img1(N,N);
    image_cuda<type> img2(N,N);
    image_cuda<type> img3(N,N);
    img1.ones();
    img2.ones();
    img3.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img3 = img1 + img2;
    };
};

static void vector_add_image_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cpu<type> img1(N,N,N);
    image_cpu<type> img2(N,N,N);
    image_cpu<type> img3(N,N,N);
    img1.ones();
    img2.ones();
    img3.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img3 = img1 + img2;
    };
};

static void vector_add_image_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_gpu<type> img1(N,N,N);
    image_gpu<type> img2(N,N,N);
    image_gpu<type> img3(N,N,N);
    img1.ones();
    img2.ones();
    img3.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img3 = img1 + img2;
    };
};

static void vector_add_image_cuda_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cuda<type> img1(N,N,N);
    image_cuda<type> img2(N,N,N);
    image_cuda<type> img3(N,N,N);
    img1.ones();
    img2.ones();
    img3.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img3 = img1 + img2;
    };
};




static void vector_assign_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_cpu<type> img1(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        img1.assign(2.3);
    };
};


// Register the function as a benchmark
BENCHMARK(vector_assign_image_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);

BENCHMARK(vector_add_image_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
// BENCHMARK(vector_add_image_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
// BENCHMARK(vector_add_image_cuda_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(vector_add_image_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
// BENCHMARK(vector_add_image_gpu_3d)->RangeMultiplier(10)->Range(10, 400);
// BENCHMARK(vector_add_image_cuda_3d)->RangeMultiplier(10)->Range(10, 400);




// Run the benchmark
BENCHMARK_MAIN();
