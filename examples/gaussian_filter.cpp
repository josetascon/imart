/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-12-05 17:46:42
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

    auto img3a = image_cpu<type>::new_pointer(5,4,4);
    img3a->random();
    auto img3b = gaussian_filter(img3a, (type)1.0);

    img3a->print_data("image 3d");
    img3b->print_data("gaussian filter to image 3d");


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

    std::string filename3 = "./examples/images/scaled_rire_ct_001.nrrd";
    auto img3i = image_cpu<short>::new_pointer(3);
    img3i->read(filename3);
    // img3i->print_data();
    
    auto img3f = image_cpu<type>::new_pointer();
    cast(*img3i, *img3f);
    // img3f->print_data();

    auto img3o = gaussian_filter(img3f, (type)1.0);
    // img3o->print_data();
    auto img3w = image_cpu<short>::new_pointer(3);
    cast(*img3o, *img3w);

    img3w->write("./scaled_rire_ct_001_blur.nrrd");


    return 0;
};