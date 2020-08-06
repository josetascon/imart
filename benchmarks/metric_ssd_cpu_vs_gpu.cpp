/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-02 14:20:33
*/

#include <iostream>
#include <benchmark/benchmark.h>

// local headers
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"
#include "../src/inearest.h"
#include "../src/ssd.h"


using namespace imart;

// *** use here linear also;

// Function to be timed
static void bm_ssd_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);

    auto image0 = image_cpu<type>::new_pointer(N,N);
    auto image1 = image_cpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_cpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();
    
    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.9, 1.3, 8.0});
    auto taffine = affine<type>::new_pointer(2, params);
    auto interp0 = inearest_cpu<type>::new_pointer(image0);
    x1 = taffine->apply(x0);
    image1 = interp0->apply(x1);

    auto tr = affine<type>::new_pointer(2);
    auto ssd1 = ssd<type,vector_cpu<type>>::new_pointer(image0, image1, tr);
    for (auto _ : state)
    {
        // This code gets timed
        ssd1->cost();
    };
};

static void bm_ssd_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_cpu<type>::new_pointer(N,N,N);
    auto image1 = image_cpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_cpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.05, 1.2, 0.03, 0, -0.04, 1, 11.324, 201.4, 8.0});
    auto taffine = affine<type>::new_pointer(3, params);
    auto interp0 = inearest_cpu<type>::new_pointer(image0);
    x1 = taffine->apply(x0);
    image1 = interp0->apply(x1);

    auto tr = affine<type>::new_pointer(3);
    auto ssd1 = ssd<type,vector_cpu<type>>::new_pointer(image0, image1, tr);
    for (auto _ : state)
    {
        // This code gets timed
        ssd1->cost();
    };
};

static void bm_ssd_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_gpu<type>::new_pointer(N,N);
    auto image1 = image_gpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_gpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_gpu<type>::pointer params(new image_gpu<type>{1.1, 0.1, -0.2, 0.9, 1.3, 8.0});
    auto taffine = affine<type,vector_ocl<type>>::new_pointer(2, params);
    auto interp0 = inearest_gpu<type>::new_pointer(image0);
    x1 = taffine->apply(x0);
    image1 = interp0->apply(x1);

    // clear data to free memory
    x0.reset();
    x1.reset();

    auto tr = affine<type,vector_ocl<type>>::new_pointer(2);
    auto ssd1 = ssd<type,vector_ocl<type>>::new_pointer(image0, image1, tr);
    for (auto _ : state)
    {
        // This code gets timed
        ssd1->cost();
    };
};

static void bm_ssd_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_gpu<type>::new_pointer(N,N,N);
    auto image1 = image_gpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_gpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_gpu<type>::pointer params(new image_gpu<type>{1.1, 0.1, -0.2, 0.05, 1.2, 0.03, 0, -0.04, 1, 11.324, 201.4, 8.0});
    auto taffine = affine<type,vector_ocl<type>>::new_pointer(3, params);
    auto interp0 = inearest_gpu<type>::new_pointer(image0);
    x1 = taffine->apply(x0);
    image1 = interp0->apply(x1);

    // clear data to free memory
    x0.reset();
    x1.reset();

    auto tr = affine<type,vector_ocl<type>>::new_pointer(3);
    auto ssd1 = ssd<type,vector_ocl<type>>::new_pointer(image0, image1, tr);
    for (auto _ : state)
    {
        // This code gets timed
        ssd1->cost();
    };
};

// Function to be timed
static void bm_ssd_derivative_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);

    auto image0 = image_cpu<type>::new_pointer(N,N);
    auto image1 = image_cpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_cpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();
    
    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.9, 1.3, 8.0});
    auto taffine = affine<type>::new_pointer(2, params);
    auto interp0 = inearest_cpu<type>::new_pointer(image0);
    x1 = taffine->apply(x0);
    image1 = interp0->apply(x1);

    auto tr = affine<type>::new_pointer(2);
    auto ssd1 = ssd<type,vector_cpu<type>>::new_pointer(image0, image1, tr);
    ssd1->cost();
    for (auto _ : state)
    {
        // This code gets timed
        ssd1->derivative();
    };
};

static void bm_ssd_derivative_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_cpu<type>::new_pointer(N,N,N);
    auto image1 = image_cpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_cpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.05, 1.2, 0.03, 0, -0.04, 1, 11.324, 201.4, 8.0});
    auto taffine = affine<type>::new_pointer(3, params);
    auto interp0 = inearest_cpu<type>::new_pointer(image0);
    x1 = taffine->apply(x0);
    image1 = interp0->apply(x1);

    auto tr = affine<type>::new_pointer(3);
    auto ssd1 = ssd<type,vector_cpu<type>>::new_pointer(image0, image1, tr);
    ssd1->cost();
    for (auto _ : state)
    {
        // This code gets timed
        ssd1->derivative();
    };
};

static void bm_ssd_derivative_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_gpu<type>::new_pointer(N,N);
    auto image1 = image_gpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_gpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_gpu<type>::pointer params(new image_gpu<type>{1.1, 0.1, -0.2, 0.9, 1.3, 8.0});
    auto taffine = affine<type,vector_ocl<type>>::new_pointer(2, params);
    auto interp0 = inearest_gpu<type>::new_pointer(image0);
    x1 = taffine->apply(x0);
    image1 = interp0->apply(x1);

    // clear data to free memory
    x0.reset();
    x1.reset();
    taffine.reset();
    interp0.reset();

    auto tr = affine<type,vector_ocl<type>>::new_pointer(2);
    auto ssd1 = ssd<type,vector_ocl<type>>::new_pointer(image0, image1, tr);
    ssd1->cost();
    for (auto _ : state)
    {
        // This code gets timed
        ssd1->derivative();
    };
};

static void bm_ssd_derivative_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_gpu<type>::new_pointer(N,N,N);
    auto image1 = image_gpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_gpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_gpu<type>::pointer params(new image_gpu<type>{1.1, 0.1, -0.2, 0.05, 1.2, 0.03, 0, -0.04, 1, 11.324, 201.4, 8.0});
    auto taffine = affine<type,vector_ocl<type>>::new_pointer(3, params);
    auto interp0 = inearest_gpu<type>::new_pointer(image0);
    x1 = taffine->apply(x0);
    image1 = interp0->apply(x1);

    // clear data to free memory
    x0.reset();
    x1.reset();
    taffine.reset();
    interp0.reset();

    auto tr = affine<type,vector_ocl<type>>::new_pointer(3);
    auto ssd1 = ssd<type,vector_ocl<type>>::new_pointer(image0, image1, tr);
    ssd1->cost();
    for (auto _ : state)
    {
        // This code gets timed
        ssd1->derivative();
    };
};

// Register the function as a benchmark
BENCHMARK(bm_ssd_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_ssd_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_ssd_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_ssd_gpu_3d)->RangeMultiplier(10)->Range(10, 400);

BENCHMARK(bm_ssd_derivative_cpu_2d)->RangeMultiplier(10)->Range(10, 5000);
BENCHMARK(bm_ssd_derivative_gpu_2d)->RangeMultiplier(10)->Range(10, 5000);
BENCHMARK(bm_ssd_derivative_cpu_3d)->RangeMultiplier(10)->Range(10, 300);
BENCHMARK(bm_ssd_derivative_gpu_3d)->RangeMultiplier(10)->Range(10, 300);


// Run the benchmark
BENCHMARK_MAIN();
