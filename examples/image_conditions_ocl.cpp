/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-29 18:25:45
*/

// File to test utilities such as: pad, unpad, normalize

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"

using namespace imart;

int main()
{
    using type = float;
    // ============================================
    //      Testing image condition function
    // ============================================
    // test condition
    auto img10 = image_ocl<type>::new_pointer(8,5);
    img10->random();
    img10->print();
    img10->print_data("random values");
    // auto img11 = *img10 < 0.5;
    img10->replace(*img10 < 0.5, 0.0);
    img10->print_data("if (img < 0.5) then equal 0");

    // auto img10 = image_cpu<char>::new_pointer(8,5);
    // img10->print();
    // img10->print_data();
    
    return 0;
};