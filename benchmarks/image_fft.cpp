/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-21 21:02:14
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
static void bm_fft_image_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::vector a = fft(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_fft_image_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cpu<type>::pointer img1 = image_cpu<type>::new_pointer(N,N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_cpu<type>::vector a = fft(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

#ifdef IMART_WITH_OPENCL
// Function to be timed
static void bm_fft_image_opencl_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::vector a = fft(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_fft_image_opencl_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_opencl<type>::pointer img1 = image_opencl<type>::new_pointer(N,N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_opencl<type>::vector a = fft(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};
#endif

#ifdef IMART_WITH_CUDA
// Function to be timed
static void bm_fft_image_cuda_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::vector a = fft(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};

// Function to be timed
static void bm_fft_image_cuda_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    image_cuda<type>::pointer img1 = image_cuda<type>::new_pointer(N,N,N);
    img1->random();
    for (auto _ : state)
    {
        // This code gets timed
        image_cuda<type>::vector a = fft(img1);
    };
    // std::cout << "Image size: " << img1.get_total_elements() << std::endl;
};
#endif

static void CustomArguments2d(benchmark::internal::Benchmark* b)
{
    std::vector<int> list = {10,20,40,80,100,200,400,800,1000,2000,4000};
    for(int i = 0; i < list.size(); i++)
        b->Args({list[i]});
};

static void CustomArguments3d(benchmark::internal::Benchmark* b)
{
    std::vector<int> list = {10,20,40,80,100,200,400};
    for(int i = 0; i < list.size(); i++)
        b->Args({list[i]});
};

// Register the function as a benchmark
// BENCHMARK(bm_fft_image_cpu_2d)->Arg(10000);
// BENCHMARK(bm_fft_image_cpu_2d)->RangeMultiplier(22)->Range(22, 484);
// BENCHMARK(bm_fft_image_cpu_3d)->RangeMultiplier(2)->Range(8, 32);
// BENCHMARK(bm_fft_image_opencl_2d)->RangeMultiplier(22)->Range(22, 484);
// BENCHMARK(bm_fft_image_opencl_3d)->RangeMultiplier(2)->Range(8, 32);
// BENCHMARK(bm_fft_image_cuda_2d)->RangeMultiplier(22)->Range(22, 484);
// BENCHMARK(bm_fft_image_cuda_3d)->RangeMultiplier(2)->Range(8, 32);

BENCHMARK(bm_fft_image_cpu_2d)->Apply(CustomArguments2d);
BENCHMARK(bm_fft_image_cpu_3d)->Apply(CustomArguments3d);

#ifdef IMART_WITH_OPENCL
BENCHMARK(bm_fft_image_opencl_2d)->Apply(CustomArguments2d);
BENCHMARK(bm_fft_image_opencl_3d)->Apply(CustomArguments3d);
#endif

#ifdef IMART_WITH_CUDA
BENCHMARK(bm_fft_image_cuda_2d)->Apply(CustomArguments2d);
BENCHMARK(bm_fft_image_cuda_3d)->Apply(CustomArguments3d);
#endif

// Run the benchmark
BENCHMARK_MAIN();
