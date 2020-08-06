/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-22 12:54:03
*/

#include <iostream>
#include <benchmark/benchmark.h>

// local headers
#include "../src/vector_vcl.h"
#include "../src/vector_ocl.h"

using namespace imart;

// Function to be timed
static void scalar_operation_vcl_vector(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    typename vector_vcl<type>::pointer input = vector_vcl<type>::new_pointer(N,2.0);
    typename vector_vcl<type>::pointer out;
    type scalar = 1.0;

    for (auto _ : state)
    {
        // This code gets timed
        int size = input->size();
        vector_vcl<type>::pointer vscale = vector_vcl<type>::new_pointer(size,scalar);
        vector_vcl<type>::pointer output = vector_vcl<type>::new_pointer(size);
        output = input->operator+(*vscale);
        out = output;
    };
    // std::cout<< out->operator[](50) << "\n"; // check result with index 50
};

static void scalar_operation_vcl_kernel(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    typename vector_vcl<type>::pointer input = vector_vcl<type>::new_pointer(N,2.0);
    typename vector_vcl<type>::pointer out;
    type scalar = 1.0;

    for (auto _ : state)
    {
        // This code gets timed
        int size = input->size();
        typename vector_vcl<type>::pointer output = vector_vcl<type>::new_pointer(size);
        std::string str_kernel = kernel_scalar( string_type<type>(), "+");
        viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
        viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
        viennacl::ocl::enqueue(kkk(*input, *output, scalar));
        out = output;
    };
    // std::cout<< out->operator[](50) << "\n"; // check result with index 50
};

static void scalar_operation_ocl_kernel(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    typename vector_ocl<type>::pointer input = vector_ocl<type>::new_pointer(N,2.0);
    typename vector_ocl<type>::pointer out;
    type scalar = 1.0;

    for (auto _ : state)
    {
        // This code gets timed
        int size = input->size();
        typename vector_ocl<type>::pointer output = vector_ocl<type>::new_pointer(size);
        std::string str_kernel = kernel_scalar( string_type<type>(), "+");
        cl_manager.program(str_kernel, "kernel_scalar");
        cl_manager.arguments(*(input->get_buffer()), *(output->get_buffer()), scalar);
        cl_manager.execute(size);
        out = output;
    };
    // std::cout<< out->operator[](50) << "\n"; // check result with index 50
};


// Register the function as a benchmark
BENCHMARK(scalar_operation_vcl_vector)->RangeMultiplier(10)->Range(10, pow(10,8));
BENCHMARK(scalar_operation_vcl_kernel)->RangeMultiplier(10)->Range(10, pow(10,8));
BENCHMARK(scalar_operation_ocl_kernel)->RangeMultiplier(10)->Range(10, pow(10,8));


// Run the benchmark
BENCHMARK_MAIN();
