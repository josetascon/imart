/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-08 06:57:28
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
    auto image2 = image_gpu<type>::new_pointer(5,3);

    // Create identity transform
    auto dfield1 = dfield_gpu<type>::new_pointer(image2);
    // Create other transform
    auto dfield2 = dfield_gpu<type>::new_pointer(image2);

    // Setup dfield2
    dfield2->get_parameters(0)->assign(2.0);
    dfield2->get_parameters(1)->assign(-0.5);

    // Print fields
    dfield1->print();
    dfield1->print_data();

    dfield2->print();
    dfield2->print_data();

    // *** TO BE IMPLEMENTED
    // // First point to transform
    // std::vector<type> point1({1.0,2.0});
    // std::vector<type> point2;

    // point2 = dfield1->apply(point1);
    // std::cout << "Input point:" << std::endl;
    // std::cout << point1[0] << " " << point1[1] << std::endl;
    // std::cout << "Transformed point (identity):" << std::endl;
    // std::cout << point2[0] << " " << point2[1] << std::endl;
    
    // // Second point to transform
    // std::vector<type> point3({1.0,2.0});
    // std::vector<type> point4;
    
    // point4 = dfield2->apply(point3);
    // std::cout << "Input point:" << std::endl;
    // std::cout << point3[0] << " " << point3[1] << std::endl;
    // std::cout << "Transformed point:" << std::endl;
    // std::cout << point4[0] << " " << point4[1] << std::endl;
    // std::cout << std::endl;
    
    // ============================================
    //      Testing 2d dfield with grids
    // ============================================
    // Grid
    auto x0 = grid_gpu<type>::new_pointer(image2);
    auto x1 = grid_gpu<type>::new_pointer();

    // Print grid
    x0->print_data("grid x0");
    
    x1 = dfield1->apply(x0);
    x1->print_data("grid x1, dfield 1 (identity)");

    x1 = dfield2->apply(x0);
    x1->print_data("grid x1, dfield 2");
    
    // ============================================
    //      Testing 3d dfield with points
    // ============================================
    // Init
    std::cout << std::string(45,'=') << std::endl;
    std::cout << std::string(15,' ') << "Example 2d" << std::endl;
    std::cout << std::string(45,'=') << std::endl;

    // Image to create grid
    auto image3 = image_gpu<type>::new_pointer(4,2,3);

    // Create identity transform
    auto dfield11 = dfield_gpu<type>::new_pointer(image3);
    // Create other transform
    auto dfield12 = dfield_gpu<type>::new_pointer(image3);

    // Setup dfield2
    dfield12->get_parameters(0)->assign(2.0);
    dfield12->get_parameters(1)->assign(-0.5);
    dfield12->get_parameters(2)->assign(1.2);

    // Print fields
    dfield11->print();
    dfield11->print_data();

    dfield12->print();
    dfield12->print_data();

    // *** TO BE IMPLEMENTED
    // // First point to transform
    // std::vector<type> point11({1.0,-2.0,8.0});
    // std::vector<type> point12;

    // point12 = dfield11->apply(point11);
    // std::cout << "Input point:" << std::endl;
    // std::cout << point11[0] << " " << point11[1] << " " << point11[2] << std::endl;
    // std::cout << "Transformed point (identity):" << std::endl;
    // std::cout << point12[0] << " " << point12[1] << " " << point12[2] << std::endl;

    // // Second point to transform
    // std::vector<type> point13({1.0,2.0,3.0});
    // std::vector<type> point14;

    // point14 = dfield12->apply(point13);
    // std::cout << "Input point:" << std::endl;
    // std::cout << point13[0] << " " << point13[1] << " " << point13[2] << std::endl;
    // std::cout << "Transformed point:" << std::endl;
    // std::cout << point14[0] << " " << point14[1] << " " << point14[2] << std::endl;
    // std::cout << std::endl;

    // ============================================
    //      Testing 3d dfield with grids
    // ============================================
    auto x10 = grid_gpu<type>::new_pointer(image3);
    auto x11 = grid_gpu<type>::new_pointer(3);
    
    x10->print_data("grid x10");

    x11 = dfield11->apply(x10);
    x11->print_data("grid x11, dfield 11 (identity)");

    x11 = dfield12->apply(x10);
    x11->print_data("grid x11, dfield 12");
    
    return 0;
};