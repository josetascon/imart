/*
* @Author: Jose Tascon
* @Date:   2019-11-07 10:13:08
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-23 18:59:23
*/


// std libs
#include <iostream>

// local libs
#include "../src/image.h"
#include "../src/image_utils.h"

using namespace imart;

int main()
{
    using type = float;

    // 2d
    std::cout<< "========================================" << std::endl;
    std::cout<< "   Two dimensional convolution with fft" << std::endl;
    std::cout<< "========================================" << std::endl;
    int w = 10;
    int h = 10;
    auto img1 = image_cpu<type>::new_pointer(w, h);
    type * p = img1->get_data()->data();
    for (int k=0; k<w*h; ++k) *(p+k) = k;
    img1->print_data("image 2d");
    
    auto kernel = image_cpu<type>::new_pointer(3,3);
    kernel->ones();
    kernel->print_data("kernel");

    auto vimg2 = convolution(img1, kernel);
    vimg2->print_data("convolution: image*kernel");
    
    auto vimg1 = fftconvolution(img1, kernel);
    vimg1->print_data("convolution fft: image*kernel");

    
    
    
    // // 3d
    // std::cout<< "========================================" << std::endl;
    // std::cout<< "          Three dimensional fft" << std::endl;
    // std::cout<< "========================================" << std::endl;
    // image_cpu<type> img11(3,4,2);
    // std::vector<type> vv(24);
    // for(int i = 0; i<vv.size(); i++) vv[i] = (type)i;
    // img11.get_data()->read_ram(vv.data(),vv.size());
    // img11.print_data("image 3d");

    // image_cpu<type>::vector vimg12;
    // vimg12 = fft(img11);
    // vimg12[0]->print_data("fft(image 3d) -> real");
    // vimg12[1]->print_data("fft(image 3d) -> img");

    // image_cpu<type>::pointer img13;
    // img13 = ifft(vimg12);
    // img13->print_data("ifft(fft(image 3d)) -> real");

    return 0;
}