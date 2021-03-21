/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-07-20 16:08:43
*/

#include <iostream>
#include <benchmark/benchmark.h>

// libs
#include "../src/vector_opencl.h"

// Testing google benchmark library
static void memory_cpu_to_opencl(benchmark::State& state)
{
    // Perform setup here
    typedef float        scalar_type;
    //typedef double    scalar_type; //use this if your GPU supports double precision
    int N = state.range(0);
    // std::cout << "Vector memory: " << 256*4 << "Bytes" << std::endl;
    std::vector<scalar_type>        std_vec(N, (scalar_type)N);
    imart::vector_opencl<scalar_type>  opencl_vec(N);
    for (auto _ : state)
    {
        // This code gets timed
        opencl_vec.read_ram(std_vec.data(), std_vec.size());
        // opencl_vec.print_data();
    };
};

// Register the function as a benchmark
// std::cout << "Vector memory: " << 256*4 << "Bytes" << std::endl;
BENCHMARK(memory_cpu_to_opencl)->RangeMultiplier(10)->Range(1000, pow(10,9)); // range loop
// BENCHMARK(memory_cpu_to_opencl)->DenseRange(1<<10,1<<30,1<<10);
// BENCHMARK(memory_cpu_to_opencl)->Arg(256);


// Run the benchmark
BENCHMARK_MAIN();
