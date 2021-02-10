/*
* @Author: Jose Tascon
* @Date:   2019-11-07 10:13:08
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-01 00:46:42
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
    std::cout<< "          Two dimensional fft" << std::endl;
    std::cout<< "========================================" << std::endl;
    image_cuda<type> img1(4,3);
    vector_cuda<type>::pointer tmp( new vector_cuda<type>({0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0}));
    img1.set_data(tmp);
    img1.print_data("image 2d");

    image_cuda<type>::vector vimg2;
    vimg2 = fft(img1);
    vimg2[0]->print_data("fft(image 2d) -> real");
    vimg2[1]->print_data("fft(image 2d) -> img");
    
    image_cuda<type>::pointer img3;
    img3 = ifft(vimg2);
    img3->print_data("ifft(fft(image 2d)) -> real");

    // 3d
    std::cout<< "========================================" << std::endl;
    std::cout<< "          Three dimensional fft" << std::endl;
    std::cout<< "========================================" << std::endl;
    image_cuda<type> img11(3,4,2);
    std::vector<type> vv(24);
    for(int i = 0; i<vv.size(); i++) vv[i] = (type)i;
    img11.get_data()->read_ram(vv.data(),vv.size());
    img11.print_data("image 3d");

    image_cuda<type>::vector vimg12;
    vimg12 = fft(img11);
    vimg12[0]->print_data("fft(image 3d) -> real");
    vimg12[1]->print_data("fft(image 3d) -> img");

    image_cuda<type>::pointer img13;
    img13 = ifft(vimg12);
    img13->print_data("ifft(fft(image 3d)) -> real");

    return 0;
}