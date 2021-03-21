/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-18 15:57:42
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
    std::string file_fixed = "./examples/images/lenag0.png";
    auto img_fixed_ = image_opencl<intype>::new_pointer();
    img_fixed_->read(file_fixed);

    std::string file_moving = "./examples/images/lenag1.png";
    auto img_moving_ = image_opencl<intype>::new_pointer();
    img_moving_->read(file_moving);

    auto img_fixed = image_opencl<type>::new_pointer();
    auto img_moving = image_opencl<type>::new_pointer();
    cast(*img_fixed_, *img_fixed);
    cast(*img_moving_, *img_moving);

    // img_fixed->print("Fixed Image");
    // std::cout << "min: " << img_fixed->min() << std::endl;
    // std::cout << "max: " << img_fixed->max() << std::endl;

    // img_moving->print("Moving Image");
    // std::cout << "min: " << img_moving->min() << std::endl;
    // std::cout << "max: " << img_moving->max() << std::endl;

    // *img_fixed = normalize<type>(*img_fixed);
    // *img_moving = normalize<type>(*img_moving);

    // img_fixed->print();
    // img_fixed->print_data();

    // img_moving->print();
    // img_moving->print_data();

    auto trfm = dfield<type,vector_opencl<type>>::new_pointer(img_fixed);
    trfm->set_sigma_elastic(2.0);
    // trfm->print();
    // trfm->print_data();

    auto demons1 = demons<type,vector_opencl<type>>::new_pointer(img_fixed, img_moving, trfm);
    auto opt = gradient_descent<type,vector_opencl<type>>::new_pointer();
    // opt->set_tolerance(1e-5);

    auto registro = registration<type,vector_opencl<type>>::new_pointer(img_fixed, img_moving, trfm);
    registro->set_metric(demons1);
    registro->set_optimizer(opt);
    // registro->set_levels_iterations(std::vector<int>{300,200,100});

    registro->apply();

    auto transformation = registro->get_transform();
    auto interpolation = ilinear_opencl<type>::new_pointer(img_moving);
    auto x0 = grid_opencl<type>::new_pointer(img_fixed);
    auto moving_warped = interpolation->apply(transformation->apply(x0));
    // moving_warped->print();
    // moving_warped->print_data();

    auto moving_cast = image_opencl<intype>::new_pointer();
    cast(*moving_warped, *moving_cast);
    moving_cast->write("./multilevel_demons_warped_opencl.png");


    return 0;
};