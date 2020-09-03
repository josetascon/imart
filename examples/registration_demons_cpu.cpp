/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-31 23:48:19
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/dfield.h"
#include "../src/demons.h"
#include "../src/gradient_descent.h"
#include "../src/registration.h"

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
    auto img_fixed_ = image_cpu<intype>::new_pointer();
    img_fixed_->read(file_fixed);

    std::string file_moving = "./examples/images/lenag2.png";
    auto img_moving_ = image_cpu<intype>::new_pointer();
    img_moving_->read(file_moving);

    auto img_fixed = image_cpu<type>::new_pointer();
    auto img_moving = image_cpu<type>::new_pointer();
    cast(*img_fixed_, *img_fixed);
    cast(*img_moving_, *img_moving);

    // img_fixed->print("Fixed Image");
    // std::cout << "min: " << img_fixed->min() << std::endl;
    // std::cout << "max: " << img_fixed->max() << std::endl;

    // img_moving->print("Moving Image");
    // std::cout << "min: " << img_moving->min() << std::endl;
    // std::cout << "max: " << img_moving->max() << std::endl;

    *img_fixed = normalize<type>(*img_fixed);
    *img_moving = normalize<type>(*img_moving);

    // img_fixed->print();
    // img_fixed->print_data();

    // img_moving->print();
    // img_moving->print_data();

    auto trfm = dfield<type,vector_cpu<type>>::new_pointer(img_moving);
    // trfm->print();
    // trfm->print_data();

    auto demons1 = demons<type,vector_cpu<type>>::new_pointer(img_fixed, img_moving, trfm);
    auto opt = gradient_descent<type,vector_cpu<type>>::new_pointer();

    auto registro = registration<type,vector_cpu<type>>::new_pointer(img_fixed, img_moving, trfm);
    registro->set_metric(demons1);
    registro->set_optimizer(opt);

    registro->apply();


    auto moving_warped = demons1->warped_moving();
    // moving_warped->print_data();

    auto moving_cast = image_cpu<intype>::new_pointer();
    cast((*moving_warped)*(type(255)), *moving_cast);

    // auto output = image_cpu<intype>::new_pointer();
    // *output = normalize<intype>(*moving_cast, min, max);

    // auto moving = (*moving_warped)*(type(255));
    // moving.print_data();
    
    // cast(*moving, *output);
    // cast(*moving_warped, *output);
    // output->print_data();
    moving_cast->write("./multilevel_demons_warped_cpu.png");
    // output->write("./demons_warped.png");


    return 0;
};