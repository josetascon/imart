/*
* @Author: Jose Tascon
* @Date:   2019-11-13 09:20:57
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2019-11-15 11:34:44
*/


// std libs
#include <iostream>

// libs
#include <benchmark/benchmark.h>

// local libs
#include "../src/image_2d.h"


// Testing google benchmark library

static void BM_Image_Init(benchmark::State& state)
{
    // Perform setup here
    
    for (auto _ : state) {
        // This code gets timed
        image_2d<double> img1(100000000,1);
        // SomeFunction();
    };
};

static void BM_Image_Zero(benchmark::State& state)
{
    // Perform setup here
    image_2d<double> img1(100000000,1);

    for (auto _ : state) {
        // This code gets timed
        img1.zeros();
        // SomeFunction();
    };
};

static void BM_Image_Random(benchmark::State& state)
{
    // Perform setup here
    image_2d<double> img1(100000000,1);

    for (auto _ : state) {
        // This code gets timed
        img1.random();
        // SomeFunction();
    };
};

static void BM_Image_Sum(benchmark::State& state)
{
    // Perform setup here
    image_2d<double> img1(100000000,2);
    image_2d<double> img2(100000000,2);
    image_2d<double> img3(100000000,2);

    for (auto _ : state) {
        // This code gets timed
        img3 = img1 + img2;
        // SomeFunction();
    };
};

// Register the function as a benchmark
// BENCHMARK(BM_Image_Init);
// BENCHMARK(BM_Image_Zero);
BENCHMARK(BM_Image_Random);
// BENCHMARK(BM_Image_Sum);

// Run the benchmark
BENCHMARK_MAIN();

// int main(int argc, char** argv)
// {
//     int num_test = 21;
//     std::vector<int> sizes({64, 128, 256, 512, 1024});


//     return 0;
// };