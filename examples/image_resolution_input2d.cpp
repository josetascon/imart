/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-07-24 10:02:18
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
    auto image1 = image_cpu<unsigned short>::new_pointer(2);
    auto image2 = image_cpu<type>::new_pointer(2);
    typename image_cpu<type>::vector vimage3(num_scales);
    typename image_cpu<unsigned short>::vector vimage4(num_scales);

    for(int k = 0; k < num_scales; k++)
    {
        vimage3[k] = image_cpu<type>::new_pointer(2);
        vimage4[k] = image_cpu<unsigned short>::new_pointer(2);
    }

    image1->read(argv[1]);
    cast(*image1,*image2);
    image2->print();
    auto mresolution_cpu = resolution<type,vector_cpu<type>>::new_pointer(image2);

    // single test
    // auto imagea = mresolution->apply(2.0);
    // auto imageo = image_cpu<unsigned short>::new_pointer(2);
    // cast(*imagea,*imageo);
    // imageo->write("./out_res_2d_cpu.png");
    
    std::string outfile1 = "./out_res_2d_cpu";
    for(int k = 0; k < num_scales; k++)
    {
        vimage3[k] = mresolution_cpu->apply(pow(2.0,k+1));
        vimage3[k]->print();
        cast(*(vimage3[k]),*(vimage4[k]));
        vimage4[k]->write(outfile1 + std::to_string(k) + ".png");
    };


    // ============================================
    //              Testing GPU
    // ============================================
    auto image11 = image_gpu<unsigned short>::new_pointer(2);
    auto image12 = image_gpu<type>::new_pointer(2);
    typename image_gpu<type>::vector vimage13(num_scales);
    typename image_gpu<unsigned short>::vector vimage14(num_scales);

    for(int k = 0; k < num_scales; k++)
    {
        vimage13[k] = image_gpu<type>::new_pointer(2);
        vimage14[k] = image_gpu<unsigned short>::new_pointer(2);
    }

    image11->read(argv[1]);
    cast(*image11,*image12);
    image12->print();
    auto mresolution_gpu = resolution<type,vector_ocl<type>>::new_pointer(image12);

    // single test
    // auto imagea1 = mresolution_gpu->apply(2.0);
    // auto imageo1 = image_gpu<unsigned short>::new_pointer(2);
    // cast(*imagea1,*imageo1);
    // imageo1->write("./out_res_2d_gpu.png");

    std::string outfile2 = "./out_res_2d_gpu";
    for(int k = 0; k < num_scales; k++)
    {
        vimage13[k] = mresolution_gpu->apply(pow(2.0,k+1));
        vimage13[k]->print();
        cast(*(vimage13[k]),*(vimage14[k]));
        vimage14[k]->write(outfile2 + std::to_string(k) + ".png");
    };

    return 0;
};