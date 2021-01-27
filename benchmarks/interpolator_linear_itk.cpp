/*
* @Author: jose
* @Date:   2019-11-13 14:27:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-24 00:21:04
*/

#include <iostream>
#include <benchmark/benchmark.h>

// local headers
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"
#include "../src/ilinear.h"

// itk
#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
// #include "itkNearestNeighborInterpolateImageFunction.h"

using namespace imart;

// *** use here linear also;

// Function to be timed
static void bm_ilinear_cpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);

    auto image0 = image_cpu<type>::new_pointer(N,N);
    auto image1 = image_cpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_cpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();
    
    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.9, 1.3, 8.0});
    auto taffine = affine<type>::new_pointer(2, params);
    auto interp0 = ilinear_cpu<type>::new_pointer(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply(x0);
        image1 = interp0->apply(x1);
        // image1 = interp0->apply_test(x1);
        
    };
    // x1.print_data();
};

static void bm_ilinear_cpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_cpu<type>::new_pointer(N,N,N);
    auto image1 = image_cpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_cpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_cpu<type>::pointer params(new image_cpu<type>{1.1, 0.1, -0.2, 0.05, 1.2, 0.03, 0, -0.04, 1, 11.324, 201.4, 8.0});
    auto taffine = affine<type>::new_pointer(3, params);
    auto interp0 = ilinear_cpu<type>::new_pointer(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply(x0);
        image1 = interp0->apply(x1);
    };
    // x1.print_data();
};

static void bm_ilinear_gpu_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_gpu<type>::new_pointer(N,N);
    auto image1 = image_gpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_gpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_gpu<type>::pointer params(new image_gpu<type>{1.1, 0.1, -0.2, 0.9, 1.3, 8.0});
    auto taffine = affine<type,vector_ocl<type>>::new_pointer(2, params);
    auto interp0 = ilinear_gpu<type>::new_pointer(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply(x0);
        image1 = interp0->apply(x1);
    };
    // x1.print_data();
};

static void bm_ilinear_gpu_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0); 

    auto image0 = image_gpu<type>::new_pointer(N,N,N);
    auto image1 = image_gpu<type>::new_pointer();
    image0->random();
    auto x0 = grid_gpu<type>::new_pointer(image0);
    auto x1 = x0->mimic();

    image_gpu<type>::pointer params(new image_gpu<type>{1.1, 0.1, -0.2, 0.05, 1.2, 0.03, 0, -0.04, 1, 11.324, 201.4, 8.0});
    auto taffine = affine<type,vector_ocl<type>>::new_pointer(3, params);
    auto interp0 = ilinear_gpu<type>::new_pointer(image0);
    
    for (auto _ : state)
    {
        // This code gets timed
        x1 = taffine->apply(x0);
        image1 = interp0->apply(x1);
    };
    // x1.print_data();
};

static void bm_ilinear_itk_2d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    const int dim = 2;
    std::vector<int> sz{N,N};

    // Image
    using ImageType = itk::Image<type, dim>;
    ImageType::Pointer image_itk = ImageType::New();
    ImageType::Pointer image_itk_result = ImageType::New();
    
    // Setup Image
    ImageType::IndexType start;
    for (int k = 0; k < sz.size(); k++) start[k] = 0;
    ImageType::SizeType sizei;
    for (int k = 0; k < sz.size(); k++) sizei[k] = sz[k];
    ImageType::RegionType region;
    region.SetSize(sizei);
    region.SetIndex(start);
    image_itk->SetRegions(region);
    image_itk->Allocate();

    // Random init
    type * p1 = image_itk->GetBufferPointer();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    for(int k=0; k<N*N; k++) *(p1+k) = (type)uniform(gen);

    // Image size to use in resample
    const ImageType::SizeType & size = image_itk->GetLargestPossibleRegion().GetSize();

    // Interpolator
    // using InterpolatorType = itk::NearestNeighborInterpolateImageFunction<ImageType, type>;
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, type>;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();

    // Transformation
    using TransformType = itk::AffineTransform<type, dim>;
    TransformType::Pointer transform = TransformType::New();

    // Set transform parameters
    TransformType::ParametersType parameters(6);
    parameters[0] = 1.1;
    parameters[1] = 0.1;
    parameters[2] = -0.2;
    parameters[3] = 0.9;
    parameters[4] = 1.3;
    parameters[5] = 8.0;
    transform->SetParameters(parameters);

    for (auto _ : state)
    {
        // This code gets timed
        // Resample
        using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType, type, type>;
        ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
        resample->SetInput(image_itk);
        resample->SetReferenceImage(image_itk);
        resample->UseReferenceImageOn();
        resample->SetSize(size);
        resample->SetDefaultPixelValue(0.0);
        resample->SetInterpolator(interpolator);
        resample->SetTransform(transform);
        resample->Update();
        image_itk_result = resample->GetOutput();  
    };
    // x1.print_data();
};

static void bm_ilinear_itk_3d(benchmark::State& state)
{
    // Perform setup here
    using type = float;     //4 Bytes
    int N = state.range(0);
    const int dim = 3;
    std::vector<int> sz{N,N,N};

    // Image
    using ImageType = itk::Image<type, dim>;
    ImageType::Pointer image_itk = ImageType::New();
    ImageType::Pointer image_itk_result = ImageType::New();
    
    // Setup Image
    ImageType::IndexType start;
    for (int k = 0; k < sz.size(); k++) start[k] = 0;
    ImageType::SizeType sizei;
    for (int k = 0; k < sz.size(); k++) sizei[k] = sz[k];
    ImageType::RegionType region;
    region.SetSize(sizei);
    region.SetIndex(start);
    image_itk->SetRegions(region);
    image_itk->Allocate();

    // Random init
    type * p1 = image_itk->GetBufferPointer();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    for(int k=0; k<N*N; k++) *(p1+k) = (type)uniform(gen);

    // Image size to use in resample
    const ImageType::SizeType & size = image_itk->GetLargestPossibleRegion().GetSize();

    // Interpolator
    // using InterpolatorType = itk::NearestNeighborInterpolateImageFunction<ImageType, type>;
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, type>;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();

    // Transformation
    using TransformType = itk::AffineTransform<type, dim>;
    TransformType::Pointer transform = TransformType::New();

    // Set transform parameters
    TransformType::ParametersType parameters(12);
    parameters[0] = 1.1;
    parameters[1] = 0.1;
    parameters[2] = -0.2;
    parameters[3] = 0.05;
    parameters[4] = 1.2;
    parameters[5] = 0.03;
    parameters[6] = 0.0;
    parameters[7] = -0.04;
    parameters[8] = 1.0;
    parameters[9] = 11.324;
    parameters[10] = 201.4;
    parameters[11] = 8.0;
    transform->SetParameters(parameters);

    for (auto _ : state)
    {
        // This code gets timed
        // Resample
        using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType, type, type>;
        ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
        resample->SetInput(image_itk);
        resample->SetReferenceImage(image_itk);
        resample->UseReferenceImageOn();
        resample->SetSize(size);
        resample->SetDefaultPixelValue(0.0);
        resample->SetInterpolator(interpolator);
        resample->SetTransform(transform);
        resample->Update();
        image_itk_result = resample->GetOutput();  
    };
    // x1.print_data();
};

// Register the function as a benchmark
// BENCHMARK(bm_ilinear_cpu_2d)->RangeMultiplier(10)->Range(10, 8000);
// BENCHMARK(bm_ilinear_gpu_2d)->RangeMultiplier(10)->Range(10, 8000);
// BENCHMARK(bm_ilinear_itk_2d)->RangeMultiplier(10)->Range(10, 8000);
// BENCHMARK(bm_ilinear_cpu_3d)->RangeMultiplier(10)->Range(10, 400);
// BENCHMARK(bm_ilinear_gpu_3d)->RangeMultiplier(10)->Range(10, 400);
// BENCHMARK(bm_ilinear_itk_3d)->RangeMultiplier(10)->Range(10, 400);

static void CustomArguments2d(benchmark::internal::Benchmark* b)
{
    // for (int i = 8; i <= 8000; i = 10*i)
    //     b->Args({i});
    std::vector<int> list = {10,20,40,80,100,200,400,800,1000,2000,4000,8000};
    for(int i = 0; i < list.size(); i++)
        b->Args({list[i]});
};

static void CustomArguments3d(benchmark::internal::Benchmark* b)
{
    // for (int i = 4; i <= 400; i = 10*i)
    //     b->Args({i});
    std::vector<int> list = {10,20,40,80,100,200,400,800};
    for(int i = 0; i < list.size(); i++)
        b->Args({list[i]});
};

BENCHMARK(bm_ilinear_cpu_2d)->Apply(CustomArguments2d);
// BENCHMARK(bm_ilinear_gpu_2d)->Apply(CustomArguments2d);
BENCHMARK(bm_ilinear_itk_2d)->Apply(CustomArguments2d);
BENCHMARK(bm_ilinear_cpu_3d)->Apply(CustomArguments3d);
// BENCHMARK(bm_ilinear_gpu_3d)->Apply(CustomArguments3d);
BENCHMARK(bm_ilinear_itk_3d)->Apply(CustomArguments3d);


// Run the benchmark
BENCHMARK_MAIN();