/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-07-20 16:13:55
*/

#include <iostream>
#include <benchmark/benchmark.h>

// libs
#include "../src/vector_opencl.h"

// Testing google benchmark library
static void memory_opencl_to_cpu(benchmark::State& state)
{
    // Perform setup here
    typedef float        scalar_type;
    //typedef double    scalar_type; //use this if your GPU supports double precision
    int N = state.range(0);
    // std::cout << "Vector memory: " << 256*4 << "Bytes" << std::endl;
    std::vector<scalar_type>        std_vec(N, (scalar_type)N);
    std::vector<scalar_type>        std_ret(N);
    imart::vector_opencl<scalar_type>  opencl_vec(N);
    opencl_vec.read_ram(std_vec.data(), std_vec.size());
    for (auto _ : state)
    {
        // This code gets timed
        opencl_vec.write_ram(std_ret.data(), opencl_vec.size());
    };
    // for(int i=0;i<std_ret.size();i++) std::cout << std_ret[i] << " ";
    // std::cout << std::endl;
};


// Register the function as a benchmark
// std::cout << "Vector memory: " << 256*4 << "Bytes" << std::endl;
BENCHMARK(memory_opencl_to_cpu)->RangeMultiplier(10)->Range(1000, pow(10,9)); // range loop
// BENCHMARK(memory_opencl_to_cpu)->DenseRange(1<<10,1<<30,1<<10);
// BENCHMARK(memory_opencl_to_cpu)->Arg(21);


// Run the benchmark
BENCHMARK_MAIN();
