/*
* @Author: Jose Tascon
* @Date:   2020-06-29 18:24:07
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-07-25 01:00:51
*/

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"
#include "../src/ilinear.h"

// itk
#include <itkImage.h>
#include <itkImageFileReader.h>
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

using namespace imart;

int main(int argc, char *argv[])
{
    using type = float;

    if( argc < 2 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " file_name.ext" << std::endl;
        return EXIT_FAILURE;
    }

    // CPU
    auto image11 = image_cpu<unsigned short>::new_pointer();
    auto image12 = image_cpu<type>::new_pointer();
    auto image13 = image_cpu<type>::new_pointer();
    auto image14 = image_cpu<unsigned short>::new_pointer();
    
    image11->read(argv[1]);
    cast(*image11, *image12);
    auto xcpu = grid_cpu<type>::new_pointer(image12);

    // image_cpu<type>::pointer params1( new image_cpu<type>({0.95, -0.1, 0.05, 0.9, 20, 0.0}) );
    image_cpu<type>::pointer params1( new image_cpu<type>({2.0, 0.0, 0.0, 1.0, 0.0, 0.0}) );
    auto taffine1 = affine<type>::new_pointer(2, params1);
    auto interp1 = ilinear<type,vector_cpu<type>>::new_pointer(image12);

    image13 = interp1->apply(taffine1->apply(xcpu));
    cast(*image13, *image14);
    image14->write("output_cpu.png");

    //GPU
    auto image21 = image_ocl<unsigned short>::new_pointer();
    auto image22 = image_ocl<type>::new_pointer();
    auto image23 = image_ocl<type>::new_pointer();
    auto image24 = image_ocl<unsigned short>::new_pointer();
    
    image21->read(argv[1]);
    cast(*image21, *image22);
    auto xgpu = grid_ocl<type>::new_pointer(image22);

    // image_ocl<type>::pointer params2( new image_ocl<type>({0.95, -0.1, 0.05, 0.9, 20, 0.0}) );
    image_ocl<type>::pointer params2( new image_ocl<type>({2.0, 0.0, 0.0, 1.0, 0.0, 0.0}) );
    auto taffine2 = affine<type,vector_ocl<type>>::new_pointer(2, params2);
    auto interp2 = ilinear<type,vector_ocl<type>>::new_pointer(image22);

    image23 = interp2->apply(taffine2->apply(xgpu));
    cast(*image23, *image24);
    image24->write("output_ocl.png");

    // ITK
    // Image
    const int dim = 2;
    using ImageType = itk::Image<unsigned short, dim>;
    using ReaderType = itk::ImageFileReader<ImageType>;
    using WriterType = itk::ImageFileWriter<ImageType>;
    
    ImageType::Pointer image_itk = ImageType::New();
    ImageType::Pointer image_out = ImageType::New();
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();


    // Read image
    reader->SetFileName(argv[1]);
    reader->Update();
    image_itk = reader->GetOutput();
    

    // Image size to use in resample
    const ImageType::SizeType & size = image_itk->GetLargestPossibleRegion().GetSize();

    // Interpolator
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();

    // Transformation
    using TransformType = itk::AffineTransform<double, dim>;
    TransformType::Pointer transform = TransformType::New();

    // Set transform parameters
    TransformType::ParametersType parameters(6);
    parameters[0] = 2.0;
    parameters[1] = 0.0;
    parameters[2] = 0.0;
    parameters[3] = 1.0;
    parameters[4] = 0.0;
    parameters[5] = 0.0;
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

    resample->Update();
    image_out = resample->GetOutput();

    // std::string out_file = 
    writer->SetFileName("./output_itk.png");
    writer->SetInput(image_out);
    writer->Update();

    return 0;
}