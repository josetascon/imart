/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-10 11:53:53
*/
\
#include <cmath>
#include <iostream>
#include <CL/cl.hpp>
#include <benchmark/benchmark.h>

// ViennaCL headers
#include <viennacl/matrix.hpp>

// Testing google benchmark library
static void memory_cpu_to_gpu(benchmark::State& state)
{
    // Perform setup here
    typedef float        scalar_type;
    //typedef double    scalar_type; //use this if your GPU supports double precision
    int N = state.range(0);
    // std::cout << "Vector memory: " << 256*4 << "Bytes" << std::endl;
    std::vector<scalar_type>      std_vec(N, (scalar_type)N);
    viennacl::vector<scalar_type> vcl_vec(N);
    
    for (auto _ : state)
    {
        // This code gets timed
        viennacl::copy(std_vec.begin(), std_vec.end(), vcl_vec.begin()); //either the STL way
    };
};


// Register the function as a benchmark
// std::cout << "Vector memory: " << 256*4 << "Bytes" << std::endl;
BENCHMARK(memory_cpu_to_gpu)->RangeMultiplier(10)->Range(1000, pow(10,9)); // range loop
// BENCHMARK(memory_cpu_to_gpu)->DenseRange(1<<10,1<<30,1<<10);
// BENCHMARK(memory_cpu_to_gpu)->Arg(256);


// Run the benchmark
BENCHMARK_MAIN();
