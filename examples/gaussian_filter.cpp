/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-08 22:17:50
*/

// File to test guassian filter

// std libs
#include <iostream>
#include <memory>
#include <complex>

// local libs
#include "../src/image.h"
#include "../src/image_utils.h"

using namespace imart;

int main()
{
    // Type
    using type = float;

    // ============================================
    //      Testing image 2d utility functions
    // ============================================
    // Kernels
    auto img2 = gaussian_kernel<type,vector_cpu<type>>(2,1,3);
    img2->print_data("Gaussian kernel 2d:");

    auto img3 = gaussian_kernel<type,vector_cpu<type>>(3,1,3);
    img3->print_data("Gaussian kernel 3d:");

    auto imga = image_cpu<type>::new_pointer(6,4);
    imga->random();
    auto imgb = gaussian_filter(imga, (type)1.0);

    imga->print_data("image");
    imgb->print_data("gaussian filter to image");


    // Input
    std::string filename = "./examples/images/cameraman20x16.tif";
    auto imgi = image_cpu<unsigned short>::new_pointer();
    imgi->read(filename);
    
    auto imgf = image_cpu<type>::new_pointer();
    cast(*imgi, *imgf);

    auto imgo = gaussian_filter(imgf, (type)1.0);
    auto imgw = image_cpu<unsigned short>::new_pointer();
    cast(*imgo, *imgw);

    imgw->write("./cameraman20x16_blur.tif");


    return 0;
};