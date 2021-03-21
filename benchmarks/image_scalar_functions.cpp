/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-18 16:23:30
*/

#include <iostream>
#include <benchmark/benchmark.h>

// local headers
#include "../src/image.h"

using namespace imart;

// =============================================================================
//                                  CPU
// =============================================================================

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

static void scalar_sub_image_cpu_2d(benchmark::State& state)
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
        img2 = img1 - scalar;
    };
};

static void scalar_mul_image_cpu_2d(benchmark::State& state)
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
        img2 = img1*scalar;
    };
};

static void scalar_div_image_cpu_2d(benchmark::State& state)
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
        img2 = img1/scalar;
    };
};

static void scalar_pow_image_cpu_2d(benchmark::State& state)
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
        img2 = img1^scalar;
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

static void scalar_sub_image_cpu_3d(benchmark::State& state)
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
        img2 = img1 - scalar;
    };
};

static void scalar_mul_image_cpu_3d(benchmark::State& state)
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
        img2 = img1*scalar;
    };
};

static void scalar_div_image_cpu_3d(benchmark::State& state)
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
        img2 = img1/scalar;
    };
};

static void scalar_pow_image_cpu_3d(benchmark::State& state)
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
        img2 = img1^scalar;
    };
};

// =============================================================================
//                                  OPENCL
// =============================================================================
#ifdef IMART_WITH_OPENCL

// Function to be timed
static void scalar_add_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N);
    image_opencl<type> img2(N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 + scalar;
    };
};

static void scalar_sub_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N);
    image_opencl<type> img2(N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 - scalar;
    };
};

static void scalar_mul_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N);
    image_opencl<type> img2(N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1*scalar;
    };
};

static void scalar_div_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N);
    image_opencl<type> img2(N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1/scalar;
    };
};

static void scalar_pow_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N);
    image_opencl<type> img2(N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1^scalar;
    };
};


static void scalar_add_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N,N);
    image_opencl<type> img2(N,N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 + scalar;
    };
};

static void scalar_sub_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N,N);
    image_opencl<type> img2(N,N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1 - scalar;
    };
};

static void scalar_mul_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N,N);
    image_opencl<type> img2(N,N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1*scalar;
    };
};

static void scalar_div_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N,N);
    image_opencl<type> img2(N,N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1/scalar;
    };
};

static void scalar_pow_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    image_opencl<type> img1(N,N,N);
    image_opencl<type> img2(N,N,N);
    type scalar = 2.0;
    img1.ones();
    img2.zeros();

    for (auto _ : state)
    {
        // This code gets timed
        img2 = img1^scalar;
    };
};
#endif

// =============================================================================
//                                  CUDA
// =============================================================================
#ifdef IMART_WITH_CUDA

// Function to be timed
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

static void scalar_sub_image_cuda_2d(benchmark::State& state)
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
        img2 = img1 - scalar;
    };
};

static void scalar_mul_image_cuda_2d(benchmark::State& state)
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
        img2 = img1*scalar;
    };
};

static void scalar_div_image_cuda_2d(benchmark::State& state)
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
        img2 = img1/scalar;
    };
};

static void scalar_pow_image_cuda_2d(benchmark::State& state)
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
        img2 = img1^scalar;
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

static void scalar_sub_image_cuda_3d(benchmark::State& state)
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
        img2 = img1 - scalar;
    };
};

static void scalar_mul_image_cuda_3d(benchmark::State& state)
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
        img2 = img1*scalar;
    };
};

static void scalar_div_image_cuda_3d(benchmark::State& state)
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
        img2 = img1/scalar;
    };
};

static void scalar_pow_image_cuda_3d(benchmark::State& state)
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
        img2 = img1^scalar;
    };
};

#endif



unsigned int end2d = 10000;
unsigned int end3d = 1000;

// Register the function as a benchmark
BENCHMARK(scalar_add_image_cpu_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_sub_image_cpu_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_mul_image_cpu_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_div_image_cpu_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_pow_image_cpu_2d)->RangeMultiplier(10)->Range(10, end2d);

#ifdef IMART_WITH_OPENCL
BENCHMARK(scalar_add_image_opencl_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_sub_image_opencl_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_mul_image_opencl_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_div_image_opencl_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_pow_image_opencl_2d)->RangeMultiplier(10)->Range(10, end2d);
#endif

#ifdef IMART_WITH_CUDA
BENCHMARK(scalar_add_image_cuda_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_sub_image_cuda_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_mul_image_cuda_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_div_image_cuda_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(scalar_pow_image_cuda_2d)->RangeMultiplier(10)->Range(10, end2d);
#endif

// 3d
BENCHMARK(scalar_add_image_cpu_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_sub_image_cpu_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_mul_image_cpu_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_div_image_cpu_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_pow_image_cpu_3d)->RangeMultiplier(10)->Range(10, end3d);

#ifdef IMART_WITH_OPENCL
BENCHMARK(scalar_add_image_opencl_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_sub_image_opencl_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_mul_image_opencl_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_div_image_opencl_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_pow_image_opencl_3d)->RangeMultiplier(10)->Range(10, end3d);
#endif

#ifdef IMART_WITH_CUDA
BENCHMARK(scalar_add_image_cuda_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_sub_image_cuda_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_mul_image_cuda_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_div_image_cuda_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(scalar_pow_image_cuda_3d)->RangeMultiplier(10)->Range(10, end3d);
#endif

// Run the benchmark
BENCHMARK_MAIN();
