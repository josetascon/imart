/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-27 06:49:35
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/dfield.h"
#include "../src/ssd.h"
#include "../src/demons.h"
#include "../src/gradient_descent.h"

using namespace imart;

int main()
{
    using intype = unsigned char;
    using type = float;
    // using type = double;
    

    // ============================================
    //          Testing metric
    // ============================================

    // Create small imame
    std::string file_fixed = "./examples/images/lenag1.png";
    auto img_fixed_ = image_opencl<intype>::new_pointer();
    img_fixed_->read(file_fixed);

    std::string file_moving = "./examples/images/lenag2.png";
    auto img_moving_ = image_opencl<intype>::new_pointer();
    img_moving_->read(file_moving);

    auto img_fixed = image_opencl<type>::new_pointer();
    auto img_moving = image_opencl<type>::new_pointer();
    cast(*img_fixed_, *img_fixed);
    cast(*img_moving_, *img_moving);

    // intype min = img_moving->min();
    // intype max = img_moving->max();

    *img_fixed = normalize<type>(*img_fixed);
    *img_moving = normalize<type>(*img_moving);

    // img_fixed->print("Fixed Image");
    // std::cout << "min: " << img_fixed->min() << std::endl;
    // std::cout << "max: " << img_fixed->max() << std::endl;

    // img_moving->print("Moving Image");
    // std::cout << "min: " << img_moving->min() << std::endl;
    // std::cout << "max: " << img_moving->max() << std::endl;

    // img_fixed->print();
    // img_fixed->print_data();

    // img_moving->print();
    // img_moving->print_data();
    
    auto trfm = dfield<type,vector_opencl<type>>::new_pointer(img_fixed);
    trfm->print();
    // trfm->print_data();
    
    auto demons1 = demons<type,vector_opencl<type>>::new_pointer(img_fixed, img_moving, trfm);
    std::cout << "cost: " << demons1->cost() << std::endl;

    auto opt = gradient_descent<type,vector_opencl<type>>::new_pointer();
    opt->set_step(1.0);
    // opt->set_iterations(1);
    opt->optimize(demons1);

    auto moving_warped = demons1->warped_moving();
    // moving_warped->print_data();

    auto moving_cast = image_opencl<intype>::new_pointer();
    cast((*moving_warped)*(type(255)), *moving_cast);

    // auto output = image_opencl<intype>::new_pointer();
    // *output = normalize<intype>(*moving_cast, min, max);

    // auto moving = (*moving_warped)*(type(255));
    // moving.print_data();
    
    // cast(*moving, *output);
    // cast(*moving_warped, *output);
    // output->print_data();
    moving_cast->write("./demons_warped_opencl.png");
    // output->write("./demons_warped.png");
    

    return 0;
};