/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-06 16:41:30
*/

// std libs
#include <iostream>
#include <memory>
#include <vector>

// local libs
#include "../src/image.h"
#include "../src/affine.h"

using namespace imart;

int main()
{
    using type = float;

    // ============================================
    //      Testing 2d affine with points
    // ============================================
    // Create identity transform
    auto affine1 = affine_opencl<type>::new_pointer(2);

    std::vector<type> point1({1.0,2.0});
    std::vector<type> point2;

    std::cout << std::string(45,'=') << std::endl;
    std::cout << std::string(15,' ') << "Example 2d" << std::endl;
    std::cout << std::string(45,'=') << std::endl;
    affine1->print();
    affine1->print_data();
    
    point2 = affine1->apply(point1);
    std::cout << "Input point:" << std::endl;
    std::cout << point1[0] << " " << point1[1] << std::endl;
    std::cout << "Transformed point (identity):" << std::endl;
    std::cout << point2[0] << " " << point2[1] << std::endl;

    // Create other transform
    image_opencl<type>::pointer params2(new image_opencl<type>{1.1, 0.5, -0.5, 0.9, 2.1, -1.1});
    auto affine2 = affine_opencl<type>::new_pointer(2,params2);

    std::vector<type> point3({1.0,2.0});
    std::vector<type> point4;

    // std::cout << affine2;
    affine2->print();
    affine2->print_data();

    point4 = affine2->apply(point3);
    std::cout << "Input point:" << std::endl;
    std::cout << point3[0] << " " << point3[1] << std::endl;
    std::cout << "Transformed point:" << std::endl;
    std::cout << point4[0] << " " << point4[1] << std::endl;
    std::cout << std::endl;
    
    // // ============================================
    // //      Testing 2d affine with grids
    // // ============================================
    auto image0 = image_opencl<type>::new_pointer(5,3);
    auto x0 = grid_opencl<type>::new_pointer(image0);
    auto x1 = grid_opencl<type>::new_pointer();
    
    x1 = affine1->apply(x0);

    x0->print_data("grid x0");
    x1->print_data("grid x1, affine 1 (identity)");

    x1 = affine2->apply(x0);
    x1->print_data("grid x1, affine 2");

    image_opencl<type>::pointer params3(new image_opencl<type>{1.0, 0.0, 0.0, 1.0, 10.0, -10.0});
    auto translation = affine_opencl<type>::new_pointer(2,params3);

    translation->print("Translation");
    translation->print_data();

    x1 = translation->apply(x0);
    x1->print_data("grid x1, translation [10, -10]");

    // ============================================
    //      Testing 3d affine with points
    // ============================================
    // Create identity transform
    auto affine11 = affine_opencl<type>::new_pointer(3);

    std::vector<type> point11({1.0,-2.0,8.0});
    std::vector<type> point12;

    std::cout << std::string(45,'=') << std::endl;
    std::cout << std::string(15,' ') << "Example 3d" << std::endl;
    std::cout << std::string(45,'=') << std::endl;
    affine11->print();
    affine11->print_data();
    
    point12 = affine11->apply(point11);
    std::cout << "Input point:" << std::endl;
    std::cout << point11[0] << " " << point11[1] << " " << point11[2] << std::endl;
    std::cout << "Transformed point (identity):" << std::endl;
    std::cout << point12[0] << " " << point12[1] << " " << point12[2] << std::endl;

    // Create other transform
    image_opencl<type>::pointer params12(new image_opencl<type>{1.1, 0.5, 0.1, -0.5, 0.9, -0.2, 0.0, 0.0, 1.0, 2.1, -1.1, 5.2});
    auto affine12 = affine_opencl<type>::new_pointer(3,params12);

    std::vector<type> point13({1.0,2.0,3.0});
    std::vector<type> point14;

    // std::cout << affine2;
    affine12->print();
    affine12->print_data();

    point14 = affine12->apply(point13);
    std::cout << "Input point:" << std::endl;
    std::cout << point13[0] << " " << point13[1] << " " << point13[2] << std::endl;
    std::cout << "Transformed point:" << std::endl;
    std::cout << point14[0] << " " << point14[1] << " " << point14[2] << std::endl;
    std::cout << std::endl;

    // ============================================
    //      Testing 3d affine with grids
    // ============================================
    auto image10 = image_opencl<type>::new_pointer(4,2,3);
    auto x10 = grid_opencl<type>::new_pointer(image10);
    auto x11 = grid_opencl<type>::new_pointer(3);
    
    x11 = affine11->apply(x10);

    x10->print_data("grid x10");
    x11->print_data("grid x11, affine 11 (identity)");

    x11 = affine12->apply(x10);
    x11->print_data("grid x11, affine 12");

    image_opencl<type>::pointer params13(new image_opencl<type>{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 10.0, -10.0 , 4.5});
    auto translation3 = affine_opencl<type>::new_pointer(3,params13);

    translation3->print("Translation");
    translation3->print_data();

    x11 = translation3->apply(x10);
    x11->print_data("grid x11, translation [10, -10, 4.5]");

    return 0;
};