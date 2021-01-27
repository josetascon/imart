/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-12-05 23:00:37
*/

// std libs
#include <iostream>
#include <memory>
#include <vector>

// local libs
#include "../src/image.h"
#include "../src/dfield.h"

using namespace imart;

int main()
{
    using type = float;

    // ============================================
    //      Testing 2d dfield with points
    // ============================================
    // Init
    std::cout << std::string(45,'=') << std::endl;
    std::cout << std::string(15,' ') << "Example 2d" << std::endl;
    std::cout << std::string(45,'=') << std::endl;

    // Image to create grid
    auto image2 = image_cpu<type>::new_pointer(5,3);

    // Create identity transform
    auto dfield2 = dfield_cpu<type>::new_pointer(image2);
    // Create other transform
    auto dfield2r = dfield_cpu<type>::new_pointer(image2);

    // Setup dfield2
    dfield2->get_parameters(0)->assign(1.5);
    dfield2->get_parameters(1)->assign(-0.2);

    // Print
    dfield2->print();
    dfield2->print_data();

    dfield2->write("./dfield2d_sample.nrrd");
    dfield2r->read("./dfield2d_sample.nrrd");

    dfield2r->print();
    dfield2r->print_data();

    // ============================================
    //      Testing 3d dfield with points
    // ============================================
    // Init
    std::cout << std::string(45,'=') << std::endl;
    std::cout << std::string(15,' ') << "Example 3d" << std::endl;
    std::cout << std::string(45,'=') << std::endl;

    // Image to create grid
    auto image3 = image_cpu<type>::new_pointer(5,3,2);

    // Create identity transform
    auto dfield3 = dfield_cpu<type>::new_pointer(image3);
    // Create other transform
    auto dfield3r = dfield_cpu<type>::new_pointer(image3);

    // Setup dfield3
    dfield3->get_parameters(0)->assign(2.0);
    dfield3->get_parameters(1)->assign(-0.5);
    dfield3->get_parameters(2)->assign(-1.0);

    // Print
    dfield3->print();
    dfield3->print_data();

    dfield3->write("./dfield3d_sample.nrrd");
    dfield3r->read("./dfield3d_sample.nrrd");

    dfield3r->print();
    dfield3r->print_data();

    return 0;
};