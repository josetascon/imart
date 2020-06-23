/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-22 12:12:23
*/
\
#include <cmath>
#include <iostream>
#include <CL/cl.hpp>
#include <benchmark/benchmark.h>

// ViennaCL headers
#include <viennacl/matrix.hpp>

// Testing google benchmark library
static void memory_gpu_to_cpu(benchmark::State& state)
{
    // Perform setup here
    typedef float        scalar_type;
    //typedef double    scalar_type; //use this if your GPU supports double precision
    int N = state.range(0);
    // std::cout << "Vector memory: " << 256*4 << "Bytes" << std::endl;
    std::vector<scalar_type>      std_vec(N, (scalar_type)N);
    std::vector<scalar_type>      std_ret(N);
    viennacl::vector<scalar_type> vcl_vec(N);
    viennacl::copy(std_vec.begin(), std_vec.end(), vcl_vec.begin());
    
    for (auto _ : state)
    {
        // This code gets timed
        viennacl::copy(vcl_vec.begin(), vcl_vec.end(), std_ret.begin());
    };
};


// Register the function as a benchmark
// std::cout << "Vector memory: " << 256*4 << "Bytes" << std::endl;
BENCHMARK(memory_gpu_to_cpu)->RangeMultiplier(10)->Range(1000, pow(10,9)); // range loop
// BENCHMARK(memory_gpu_to_cpu)->DenseRange(1<<10,1<<30,1<<10);
// BENCHMARK(memory_gpu_to_cpu)->Arg(10);


// Run the benchmark
BENCHMARK_MAIN();
