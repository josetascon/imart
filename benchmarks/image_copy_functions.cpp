/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-07-20 15:44:58
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

// Conclusion: CPU time is spent in allocation.

// Function to be timed
// static void bm_create_vector_2d(benchmark::State& state)
// {
//     // Perform setup here
//     using type = float;     //4 Bytes
//     int N = state.range(0);
//     for (auto _ : state)
//     {
//         // This code gets timed
//         std::shared_ptr<std::vector<type>> data = std::make_shared<std::vector<type>>(N*N);;
//     };
//     // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
// };

// Function to be timed
static void bm_create_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_create_image_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N,N);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_create_image_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_create_image_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N,N);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_clone_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::pointer img2 = img1->clone();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_clone_image_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::pointer img2 = img1->clone();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_clone_image_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::pointer img2 = img1->clone();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_clone_image_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::pointer img2 = img1->clone();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_copy_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::pointer img2 = img1->copy();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_copy_image_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::pointer img2 = img1->copy();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_copy_image_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::pointer img2 = img1->copy();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_copy_image_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::pointer img2 = img1->copy();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_mimic_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::pointer img2 = img1->mimic();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_mimic_image_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::pointer img2 = img1->mimic();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_mimic_image_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::pointer img2 = img1->mimic();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_mimic_image_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_gpu<type>::pointer img1 = image_gpu<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_gpu<type>::pointer img2 = img1->mimic();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};


// Register the function as a benchmark
// BENCHMARK(bm_create_vector_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_create_image_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_create_image_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_create_image_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_create_image_gpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_clone_image_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_clone_image_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_clone_image_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_clone_image_gpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_copy_image_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_copy_image_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_copy_image_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_copy_image_gpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_mimic_image_cpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_mimic_image_gpu_2d)->RangeMultiplier(10)->Range(10, 10000);
BENCHMARK(bm_mimic_image_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
BENCHMARK(bm_mimic_image_gpu_3d)->RangeMultiplier(10)->Range(10, 400);


// Run the benchmark
BENCHMARK_MAIN();
