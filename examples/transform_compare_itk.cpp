/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-03-05 18:19:50
*/


// std libs
#include <iostream>
#include <memory>

// itk
#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"
#include "../src/interpolate.h"
#include "../src/utils/timer.h"

using type = double;
const int dim = 2;

itk::Image<type,dim>::Pointer create_image(int w, int h)
{
    using ImageType = itk::Image<type, dim>;
    ImageType::Pointer image = ImageType::New();

    ImageType::IndexType start;
    start[0] = 0; // first index on X
    start[1] = 0; // first index on Y

    ImageType::SizeType size;
    size[0] = w; // size along X
    size[1] = h; // size along Y

    ImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    image->SetRegions(region);
    image->Allocate();
    return image;
};

void copy_image(itk::Image<type,dim>::Pointer input1, image<type> & input2)
{
    int num = input2.get_total_elements();
    type * p1 = input1->GetBufferPointer();
    type * p2 = input2.ptr();

    for(int k=0; k<num; k++)
    {
         *(p1+k) = *(p2+k);
    };
};


// Results

// Multicore OMP 8 threads
// Time project: 0.566975s
// Time itk: 0.188898s

// Single Core
// Time project: 0.273303s = 0.035s transform + 0.235s interp!!
// Time itk: 0.100865s


int main()
{
    // ============================================
    //          Testing resample
    // ============================================

    // Create images
    int w = 800;
    int h = 600;
    image<type> image0(w,h);
    image0.random();
    
    // Create grid
    grid<type> x0(image0);

    // Create transform
    image<type>::pointer params( new image<type>({0.9, 0.1, 0.1, 0.9, 100.5, -50.0}) );
    affine<type> affine1(2, params);

    // Create interpolation
    interpolate<type> image0_itp(image0, x0);

    // Transform
    timer t1;
    t1.start();
    // grid<type> x1 = affine1*x0;
    // image0_itp*x1;
    // image<type> image1 = image0_itp*x1;
    image<type> image1 = image0_itp*(affine1*x0);
    t1.finish();

    // x1.print_data();
    image0.print("input image");
    image1.print("transformed image");
    std::cout << std::endl;

    std::cout << "Time project: " << t1.get_time() << t1.get_units() << std::endl;

    // Image
    using ImageType = itk::Image<type, dim>;
    ImageType::Pointer image_itk;
    ImageType::Pointer image_itk_result;
    image_itk = create_image(w, h);

    copy_image(image_itk, image0);

    // Image size to use in resample
    const ImageType::SizeType & size = image_itk->GetLargestPossibleRegion().GetSize();

    // Interpolator
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, type>;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();

    // Transformation
    using TransformType = itk::AffineTransform<type, dim>;
    TransformType::Pointer transform = TransformType::New();

    // Set transform parameters
    TransformType::ParametersType parameters(6);
    parameters[0] = 0.9;
    parameters[1] = 0.1;
    parameters[2] = 0.1;
    parameters[3] = 0.9;
    parameters[4] = 100.5;
    parameters[5] = -50.0;
    transform->SetParameters(parameters);

    // Resample
    using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
    ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
    resample->SetInput(image_itk);
    resample->SetReferenceImage(image_itk);
    resample->UseReferenceImageOn();
    resample->SetSize(size);
    resample->SetDefaultPixelValue(0.0);
    resample->SetInterpolator(interpolator);
    resample->SetTransform(transform);

    // Run
    t1.start();
    resample->Update();
    resample->GetOutput();
    t1.finish();

    std::cout << "Time itk: " << t1.get_time() << t1.get_units() << std::endl;

    // ## CHECK IMAGE OUTPUT in FILES Enero 29 2020

    return 0;

};