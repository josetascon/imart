/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-10 16:16:21
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/resolution.h"

using namespace imart;

int main(int argc, char *argv[])
{
    // ============================================
    //              Testing resolution
    // ============================================
    using type = float;
    // using type = double;

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
    auto image1 = image_cpu<short>::new_pointer(3);
    auto image2 = image_cpu<type>::new_pointer(3);
    typename image_cpu<type>::vector vimage3(num_scales);
    typename image_cpu<short>::vector vimage4(num_scales);

    for(int k = 0; k < num_scales; k++)
    {
        vimage3[k] = image_cpu<type>::new_pointer(3);
        vimage4[k] = image_cpu<short>::new_pointer(3);
    }

    image1->read(argv[1]);
    cast(*image1,*image2);
    image2->print();
    auto mresolution_cpu = resolution<type,vector_cpu<type>>::new_pointer(image2);

    // single test
    // auto imagea = mresolution->apply(2.0);
    // auto imageo = image_cpu<short>::new_pointer(2);
    // cast(*imagea,*imageo);
    // imageo->write("./out_res_3d_cpu.nrrd");
    
    std::string outfile1 = "./out_res_3d_cpu";
    for(int k = 0; k < num_scales; k++)
    {
        vimage3[k] = mresolution_cpu->apply(pow(2.0,k+1));
        vimage3[k]->print();
        cast(*(vimage3[k]),*(vimage4[k]));
        vimage4[k]->write(outfile1 + std::to_string(k) + ".nrrd");
    };


    // ============================================
    //              Testing GPU
    // ============================================
    auto image11 = image_ocl<short>::new_pointer(3);
    auto image12 = image_ocl<type>::new_pointer(3);
    typename image_ocl<type>::vector vimage13(num_scales);
    typename image_ocl<short>::vector vimage14(num_scales);

    for(int k = 0; k < num_scales; k++)
    {
        vimage13[k] = image_ocl<type>::new_pointer(3);
        vimage14[k] = image_ocl<short>::new_pointer(3);
    }

    image11->read(argv[1]);
    cast(*image11,*image12);
    // image12->print();
    // image12->write("./out_res_3d_ocl.nrrd");
    auto mresolution_ocl = resolution<type,vector_ocl<type>>::new_pointer(image12);

    // single test
    // auto imagea1 = mresolution_ocl->apply(2.0);
    // auto imageo1 = image_ocl<short>::new_pointer(3);
    // cast(*imagea1,*imageo1);
    // imageo1->write("./out_res_3d_ocl.nrrd");

    std::string outfile2 = "./out_res_3d_ocl";
    for(int k = 0; k < num_scales; k++)
    {
        vimage13[k] = mresolution_ocl->apply(pow(2.0,k+1));
        vimage13[k]->print();
        cast(*(vimage13[k]),*(vimage14[k]));
        vimage14[k]->write(outfile2 + std::to_string(k) + ".nrrd");
    };

    return 0;
};