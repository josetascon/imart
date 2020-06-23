/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-22 13:00:25
*/

#include <iostream>
#include <CL/cl.hpp>
#include <benchmark/benchmark.h>

// viennacl headers
#include <viennacl/matrix.hpp>

// local headers
#include "../src/vector_vcl.h"
#include "../src/vector_ocl.h"

// Function to be timed
static void bm_create_vector_lib(benchmark::State& state)
{
    // Perform setup here
    typedef float        scalar_type; //4 Bytes
    int N = state.range(0);
    // std::cout<< "Vector memory: " << N*4 << "Bytes" << std::endl; 
    // viennacl::vector<scalar_type> vv(N);
    for (auto _ : state)
    {
        // This code gets timed
        viennacl::vector<scalar_type> vec1(N);
    };
    // std::cout << "Vector size: " << vv.size() << std::endl;
};

// Function to be timed
static void bm_create_vector_vcl(benchmark::State& state)
{
    // Perform setup here
    typedef float        scalar_type; //4 Bytes
    int N = state.range(0);
    // std::cout<< "Vector memory: " << N*4 << "Bytes" << std::endl; 
    // imart::vector_vcl<scalar_type> vv(N);
    for (auto _ : state)
    {
        // This code gets timed
        imart::vector_vcl<scalar_type> vec1(N);
    };
    // std::cout << "Vector size: " << vv.size() << std::endl;
};

// Function to be timed
static void bm_create_vector_ocl(benchmark::State& state)
{
    // Perform setup here
    typedef float        scalar_type; //4 Bytes
    int N = state.range(0);
    // std::cout<< "Vector memory: " << N*4 << "Bytes" << std::endl; 
    // imart::vector_vcl<scalar_type> vv(N);
    for (auto _ : state)
    {
        // This code gets timed
        imart::vector_ocl<scalar_type> vec1(N);
    };
    // std::cout << "Vector size: " << vv.size() << std::endl;
};

// Function to be timed
static void bm_mimic_vector_vcl(benchmark::State& state)
{
    // Perform setup here
    typedef float        scalar_type; //4 Bytes
    int N = state.range(0);
    // std::cout<< "Vector memory: " << N*4 << "Bytes" << std::endl; 
    imart::vector_vcl<scalar_type> vv(N);
    for (auto _ : state)
    {
        // This code gets timed
        // imart::vector_vcl<scalar_type>::pointer vec1 = imart::vector_vcl<scalar_type>::new_pointer(vv.size());
        imart::vector_vcl<scalar_type>::pointer vec1 = vv.mimic();
    };
    // std::cout << "Vector size: " << vv.size() << std::endl;
};

// Function to be timed
static void bm_mimic_vector_ocl(benchmark::State& state)
{
    // Perform setup here
    typedef float        scalar_type; //4 Bytes
    int N = state.range(0);
    // std::cout<< "Vector memory: " << N*4 << "Bytes" << std::endl; 
    imart::vector_ocl<scalar_type> vv(N);
    for (auto _ : state)
    {
        // This code gets timed
        // imart::vector_vcl<scalar_type>::pointer vec1 = imart::vector_vcl<scalar_type>::new_pointer(vv.size());
        imart::vector_ocl<scalar_type>::pointer vec1 = vv.mimic();
    };
    // std::cout << "Vector size: " << vv.size() << std::endl;
};

int input = 1<<18;

// Register the function as a benchmark
BENCHMARK(bm_create_vector_lib)->Arg(input);
BENCHMARK(bm_create_vector_vcl)->Arg(input);
BENCHMARK(bm_create_vector_ocl)->Arg(input);
BENCHMARK(bm_mimic_vector_vcl)->Arg(input);
BENCHMARK(bm_mimic_vector_ocl)->Arg(input);


// Run the benchmark
BENCHMARK_MAIN();
