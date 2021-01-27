/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-21 14:50:22
*/

// std libs
#include <iostream>
#include <memory>
#include <vector>

// benhcmark header
#include <benchmark/benchmark.h>

// local headers
#include "../src/image.h"
#include "../src/image_utils.h"

using namespace imart;

// Function to be timed
static void bm_gradient_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::vector a = gradient(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_gradient_image_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::vector a = gradient(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_gradient_image_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::vector a = gradient(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_gradient_image_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::vector a = gradient(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Register the function as a benchmark
// BENCHMARK(bm_gradient_image_cpu_2d)->Arg(10000);
BENCHMARK(bm_gradient_image_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_gradient_image_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_gradient_image_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_gradient_image_gpu_3d)->RangeMultiplier(10)->Range(10, 400);

// Run the benchmark
BENCHMARK_MAIN();
