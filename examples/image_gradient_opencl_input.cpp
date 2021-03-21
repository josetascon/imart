/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-21 14:49:17
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/image_utils.h"

using namespace imart;

int main()
{
    using type = float;
    // using type = double;

    // ============================================
    //          Testing metric
    // ============================================

    // Create small imame
    std::string filename = "./examples/images/cameraman20x16.tif";
    auto img = image_opencl<unsigned short>::new_pointer();
    img->read(filename);
    
    auto image0 = image_opencl<type>::new_pointer();
    cast(*img, *image0);
    image0->print_data("Image Input");

    auto v = fft(image0);
    // v[0]->print_data();

    int d = image0->get_dimension();
    std::vector<int> none(d);
    std::vector<int> extra(d);
    for(int i = 0; i < d; i++)
    {
        none[i] = 0;
        extra[i] = 2;
    };

    // ref
    // image_opencl<type> image0_pad = pad(*image0, none, extra);
    // image0->print();
    // image0_pad.print();

    // auto w = fft(image0_pad);
    // w[0]->print_data();

    //pointer
    image_opencl<type>::pointer image0_pad = pad(image0, none, extra);
    image0->print();
    image0_pad->print();

    auto w = fft(image0_pad);


    typename image_opencl<type>::vector grad(image0->get_dimension());
    grad = gradient(image0);
    // image0->print_data("i");
    grad[0]->print();
    grad[0]->print_data("di/dx");
    // grad[1]->print_data("di/dy");

    
    return 0;
};