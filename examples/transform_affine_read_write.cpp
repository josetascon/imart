/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-08 09:55:41
*/

// std libs
#include <iostream>
#include <memory>
#include <vector>

// local libs
#include "../src/affine.h"

using namespace imart;

int main()
{
    using type = float;

    // Create transform and store in file
    auto affine1 = affine_cpu<type>::new_pointer(2);
    std::initializer_list<type> list{1.2, 0.1, -0.2, 1.0, 53.1, -38.4};
    auto params = image_cpu<type>::new_pointer(list);
    affine1->set_parameters(params);
    affine1->write("./affine2d.mat");

    // Read transform back
    auto affine2 = affine_cpu<type>::new_pointer(2);
    affine2->read("./affine2d.mat");

    // Print info
    affine1->print();
    affine1->print_data();

    affine2->print();
    affine2->print_data();
    
    auto affine3 = affine_cpu<type>::new_pointer(3);
    std::initializer_list<type> list3{1.2, 0.1, -0.2, 0.0, 1.0, 0.05, 0.15, 0.03, 0.95, 53.1, -38.4, 14.32};
    auto params3 = image_cpu<type>::new_pointer(list3);
    affine3->set_parameters(params3);
    affine3->write("./affine3d.mat");

    // Read transform back
    auto affine4 = affine_cpu<type>::new_pointer(3);
    affine4->read("./affine3d.mat");

    // Print info
    affine3->print();
    affine3->print_data();

    affine4->print();
    affine4->print_data();

    return 0;
};