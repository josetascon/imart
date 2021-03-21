/*
* @Author: Jose Tascon
* @Date:   2020-06-29 18:24:07
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-08 23:45:35
*/

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/dfield.h"
#include "../src/ilinear.h"

// itk
// #include <itkImage.h>
// #include <itkImageFileReader.h>
// #include "itkAffineTransform.h"
// #include "itkResampleImageFilter.h"
// #include "itkLinearInterpolateImageFunction.h"

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
    
    auto tdfield1 = dfield_cpu<type>::new_pointer(image12);
    tdfield1->get_parameters(0)->assign(15);
    auto interp1 = ilinear<type,vector_cpu<type>>::new_pointer(image12);

    image13 = interp1->apply(tdfield1->apply(xcpu));
    cast(*image13, *image14);
    image14->write("output_dfield_cpu.png");

    //GPU
    auto image21 = image_opencl<unsigned short>::new_pointer();
    auto image22 = image_opencl<type>::new_pointer();
    auto image23 = image_opencl<type>::new_pointer();
    auto image24 = image_opencl<unsigned short>::new_pointer();
    
    image21->read(argv[1]);
    cast(*image21, *image22);
    auto xgpu = grid_opencl<type>::new_pointer(image22);

    auto tdfield2 = dfield_opencl<type>::new_pointer(image22);
    tdfield2->get_parameters(0)->assign(15);
    auto interp2 = ilinear<type,vector_opencl<type>>::new_pointer(image22);

    image23 = interp2->apply(tdfield2->apply(xgpu));
    cast(*image23, *image24);
    image24->write("output_dfield_opencl.png");

    return 0;
}