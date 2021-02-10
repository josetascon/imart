/*
* @Author: Jose Tascon
* @Date:   2020-06-29 18:24:07
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-11 00:37:06
*/

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/resolution.h"
#include "../src/icubic.h"

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

    int num_scales = 3;

    // ============================================
    //              Testing CPU
    // ============================================
    auto image1 = image_cpu<unsigned char>::new_pointer(2);
    auto image2 = image_cpu<type>::new_pointer(2);
    typename image_cpu<type>::vector vimage3(num_scales);
    typename image_cpu<unsigned char>::vector vimage4(num_scales);

    for(int k = 0; k < num_scales; k++)
    {
        vimage3[k] = image_cpu<type>::new_pointer(2);
        vimage4[k] = image_cpu<unsigned char>::new_pointer(2);
    }

    image1->read(argv[1]);
    cast(*image1,*image2);
    image2->print();

    auto mresolution_cpu = resolution<type,vector_cpu<type>>::new_pointer(image2);
    // auto interpolate = ilinear<type,vector_cpu<type>>::new_pointer(2);
    auto interpolate = icubic<type,vector_cpu<type>>::new_pointer(2);
    mresolution_cpu->set_interpolator(interpolate);

    // single test
    // auto imagea = mresolution->apply(2.0);
    // auto imageo = image_cpu<unsigned char>::new_pointer(2);
    // cast(*imagea,*imageo);
    // imageo->write("./out_res_2d_cpu.png");
    
    std::string outfile1 = "./scaleup_icubic_cpu";
    for(int k = 0; k < num_scales; k++)
    {
        vimage3[k] = mresolution_cpu->apply(pow(2.0,-1.0*(k+1)));
        vimage3[k]->print();
        cast(*(vimage3[k]),*(vimage4[k]));
        vimage4[k]->write(outfile1 + std::to_string(k) + ".png");
    };

    // ============================================
    //              Testing GPU
    // ============================================
    auto image11 = image_ocl<unsigned char>::new_pointer(2);
    auto image12 = image_ocl<type>::new_pointer(2);
    typename image_ocl<type>::vector vimage13(num_scales);
    typename image_ocl<unsigned char>::vector vimage14(num_scales);

    for(int k = 0; k < num_scales; k++)
    {
        vimage13[k] = image_ocl<type>::new_pointer(2);
        vimage14[k] = image_ocl<unsigned char>::new_pointer(2);
    }

    image11->read(argv[1]);
    cast(*image11,*image12);
    image12->print();

    auto mresolution_cpu1 = resolution<type,vector_ocl<type>>::new_pointer(image12);
    // auto interpolate1 = ilinear<type,vector_ocl<type>>::new_pointer(2);
    auto interpolate1 = icubic<type,vector_ocl<type>>::new_pointer(2);
    mresolution_cpu1->set_interpolator(interpolate1);

    // single test
    // auto imagea = mresolution->apply(2.0);
    // auto imageo = image_ocl<unsigned char>::new_pointer(2);
    // cast(*imagea,*imageo);
    // imageo->write("./out_res_2d_cpu.png");
    
    std::string outfile11 = "./scaleup_icubic_ocl";
    for(int k = 0; k < num_scales; k++)
    {
        vimage13[k] = mresolution_cpu1->apply(pow(2.0,-1.0*(k+1)));
        vimage13[k]->print();
        cast(*(vimage13[k]),*(vimage14[k]));
        vimage14[k]->write(outfile11 + std::to_string(k) + ".png");
    };

    


    // ITK
    // Image
    // const int dim = 2;
    // using ImageType = itk::Image<unsigned short, dim>;
    // using ReaderType = itk::ImageFileReader<ImageType>;
    // using WriterType = itk::ImageFileWriter<ImageType>;
    
    // ImageType::Pointer image_itk = ImageType::New();
    // ImageType::Pointer image_out = ImageType::New();
    // ReaderType::Pointer reader = ReaderType::New();
    // WriterType::Pointer writer = WriterType::New();


    // // Read image
    // reader->SetFileName(argv[1]);
    // reader->Update();
    // image_itk = reader->GetOutput();
    

    // // Image size to use in resample
    // const ImageType::SizeType & size = image_itk->GetLargestPossibleRegion().GetSize();

    // // Interpolator
    // using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
    // InterpolatorType::Pointer interpolator = InterpolatorType::New();

    // // Transformation
    // using TransformType = itk::AffineTransform<double, dim>;
    // TransformType::Pointer transform = TransformType::New();

    // // Set transform parameters
    // TransformType::ParametersType parameters(6);
    // parameters[0] = 2.0;
    // parameters[1] = 0.0;
    // parameters[2] = 0.0;
    // parameters[3] = 1.0;
    // parameters[4] = 0.0;
    // parameters[5] = 0.0;
    // transform->SetParameters(parameters);

    // // Resample
    // using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
    // ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
    // resample->SetInput(image_itk);
    // resample->SetReferenceImage(image_itk);
    // resample->UseReferenceImageOn();
    // resample->SetSize(size);
    // resample->SetDefaultPixelValue(0.0);
    // resample->SetInterpolator(interpolator);
    // resample->SetTransform(transform);

    // resample->Update();
    // image_out = resample->GetOutput();

    // // std::string out_file = 
    // writer->SetFileName("./output_cubic_itk.png");
    // writer->SetInput(image_out);
    // writer->Update();

    return 0;
}