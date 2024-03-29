/*
* @Author: Jose Tascon
* @Date:   2019-11-13 09:20:57
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-01-24 16:20:30
*/


// std libs
#include <iostream>

// libs
#include <benchmark/benchmark.h>

// local libs
#include "../src/image.h"


// Testing google benchmark library

static void BM_Image_Init(benchmark::State& state)
{
    // Perform setup here
    
    for (auto _ : state) {
        // This code gets timed
        image<double> img1(100000000,1);
        // SomeFunction();
    };
};

static void BM_Image_Zero(benchmark::State& state)
{
    // Perform setup here
    image<double> img1(100000000,1);

    for (auto _ : state) {
        // This code gets timed
        img1.zeros();
        // SomeFunction();
    };
};

static void BM_Image_Random(benchmark::State& state)
{
    // Perform setup here
    image<double> img1(100000000,1);

    for (auto _ : state) {
        // This code gets timed
        img1.random();
        // SomeFunction();
    };
};

static void BM_Image_Equal(benchmark::State& state)
{
    // Perform setup here
    image<double> img1(100000000,1);
    image<double> img2(100000000,1);
    img1.ones();

    for (auto _ : state) {
        // This code gets timed
        img2 = img1;
        // SomeFunction();
    };
};

static void BM_Image_Sum(benchmark::State& state)
{
    // Perform setup here
    image<double> img1(100000000,2);
    image<double> img2(100000000,2);
    image<double> img3(100000000,2);

    for (auto _ : state) {
        // This code gets timed
        img3 = img1 + img2;
        // SomeFunction();
    };
};


static void BM_Image_Matrix_Mult(benchmark::State& state)
{
    // Perform setup here
    image<double> img1(2000,1000);
    image<double> img2(1000,2000);
    image<double> img3(1000,1000);

    for (auto _ : state) {
        // This code gets timed
        img3 = img1._x_(img2);
        // SomeFunction();
    };
};

// Register the function as a benchmark
// BENCHMARK(BM_Image_Init);
BENCHMARK(BM_Image_Zero);
// BENCHMARK(BM_Image_Random);
// BENCHMARK(BM_Image_Equal);
// BENCHMARK(BM_Image_Sum);
// BENCHMARK(BM_Image_Matrix_Mult);

// Run the benchmark
BENCHMARK_MAIN();

// int main(int argc, char** argv)
// {
//     int num_test = 21;
//     std::vector<int> sizes({64, 128, 256, 512, 1024});


//     return 0;
// };