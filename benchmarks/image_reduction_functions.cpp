/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-18 13:21:55
*/

// std libs
#include <iostream>
#include <memory>
#include <vector>

// benhcmark header
#include <benchmark/benchmark.h>

// local headers
#include "../src/image.h"

using namespace imart;

// Function to be timed
static void bm_sum_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N);
    img1->ones();
    for (auto _ : state)
    {
        // This code gets timed
        type a = img1->sum();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_sum_image_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N,N);
    img1->ones();
    for (auto _ : state)
    {
        // This code gets timed
        type a = img1->sum();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_sum_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N);
    img1->ones();
    for (auto _ : state)
    {
        // This code gets timed
        type a = img1->sum();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_sum_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N,N);
    img1->ones();
    for (auto _ : state)
    {
        // This code gets timed
        type a = img1->sum();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Register the function as a benchmark
BENCHMARK(bm_sum_image_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_sum_image_opencl_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_sum_image_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_sum_image_opencl_3d)->RangeMultiplier(10)->Range(10, 400);

// Run the benchmark
BENCHMARK_MAIN();
