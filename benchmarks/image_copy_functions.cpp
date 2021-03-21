/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-18 15:53:31
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


// =============================================================================
//                                  CPU
// =============================================================================

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


// =============================================================================
//                                  OPENCL
// =============================================================================
#ifdef IMART_WITH_OPENCL

// Function to be timed
static void bm_create_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_create_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N,N);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_clone_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::pointer img2 = img1->clone();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_clone_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::pointer img2 = img1->clone();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_copy_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::pointer img2 = img1->copy();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_copy_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::pointer img2 = img1->copy();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_mimic_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::pointer img2 = img1->mimic();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_mimic_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::pointer img2 = img1->mimic();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};
#endif

// =============================================================================
//                                  CUDA
// =============================================================================
#ifdef IMART_WITH_CUDA

// Function to be timed
static void bm_create_image_cuda_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_create_image_cuda_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N,N);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_clone_image_cuda_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::pointer img2 = img1->clone();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_clone_image_cuda_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::pointer img2 = img1->clone();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_copy_image_cuda_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::pointer img2 = img1->copy();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_copy_image_cuda_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::pointer img2 = img1->copy();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_mimic_image_cuda_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::pointer img2 = img1->mimic();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_mimic_image_cuda_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N,N);
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::pointer img2 = img1->mimic();
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};
#endif




unsigned int end2d = 10000;
unsigned int end3d = 1000;

// Register the function as a benchmark
BENCHMARK(bm_create_image_cpu_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(bm_clone_image_cpu_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(bm_copy_image_cpu_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(bm_mimic_image_cpu_2d)->RangeMultiplier(10)->Range(10, end2d);

#ifdef IMART_WITH_OPENCL
BENCHMARK(bm_create_image_opencl_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(bm_clone_image_opencl_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(bm_copy_image_opencl_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(bm_mimic_image_opencl_2d)->RangeMultiplier(10)->Range(10, end2d);
#endif

#ifdef IMART_WITH_CUDA
BENCHMARK(bm_create_image_cuda_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(bm_clone_image_cuda_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(bm_copy_image_cuda_2d)->RangeMultiplier(10)->Range(10, end2d);
BENCHMARK(bm_mimic_image_cuda_2d)->RangeMultiplier(10)->Range(10, end2d);
#endif

BENCHMARK(bm_create_image_cpu_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(bm_clone_image_cpu_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(bm_copy_image_cpu_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(bm_mimic_image_cpu_3d)->RangeMultiplier(10)->Range(10, end3d);

#ifdef IMART_WITH_OPENCL
BENCHMARK(bm_create_image_opencl_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(bm_clone_image_opencl_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(bm_copy_image_opencl_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(bm_mimic_image_opencl_3d)->RangeMultiplier(10)->Range(10, end3d);
#endif

#ifdef IMART_WITH_CUDA
BENCHMARK(bm_create_image_cuda_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(bm_clone_image_cuda_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(bm_copy_image_cuda_3d)->RangeMultiplier(10)->Range(10, end3d);
BENCHMARK(bm_mimic_image_cuda_3d)->RangeMultiplier(10)->Range(10, end3d);
#endif

// Run the benchmark
BENCHMARK_MAIN();
