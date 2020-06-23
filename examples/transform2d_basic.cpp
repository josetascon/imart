/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-18 12:46:38
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
    // ============================================
    //      Testing affine with points
    // ============================================
    // Create identity transform
    // auto params1 = image<float>::new_pointer({1.0, 0.0, 0.0, 1.0, 0.0, 0.0});
    image<float>::pointer params1(new image<float>{1.0, 0.0, 0.0, 1.0, 0.0, 0.0});
    affine<float> affine1(2,params1);

    std::vector<float> point1({1.0,2.0});
    std::vector<float> point2;

    // std::cout << affine1;
    affine1.print();
    affine1.print_data();
    
    point2 = affine1.apply(point1);
    std::cout << "Input point:" << std::endl;
    std::cout << point1[0] << " " << point1[1] << std::endl;
    std::cout << "Transformed point (identity):" << std::endl;
    std::cout << point2[0] << " " << point2[1] << std::endl;

    // Create other transform
    image<float>::pointer params2(new image<float>{1.1, 0.5, -0.5, 0.9, 2.1, -1.1});
    affine<float> affine2(2,params2);

    std::vector<float> point3({1.0,2.0});
    std::vector<float> point4;

    // std::cout << affine2;
    affine2.print();
    affine2.print_data();

    point4 = affine2.apply(point3);
    std::cout << "Input point:" << std::endl;
    std::cout << point3[0] << " " << point3[1] << std::endl;
    std::cout << "Transformed point:" << std::endl;
    std::cout << point4[0] << " " << point4[1] << std::endl;
    std::cout << std::endl;
    
    // // ============================================
    // //      Testing affine with grids
    // // ============================================
    image<float> image0(5,3);
    grid<float> x0(image0);
    grid<float> x1;
    
    x1 = affine1.apply(x0);

    x0.print_data("grid x0");
    x1.print_data("grid x1, affine 1 (identity)");

    x1 = affine2.apply(x0);
    x1.print_data("grid x1, affine 2");

    image<float>::pointer params3(new image<float>{1.0, 0.0, 0.0, 1.0, 10.0, -10.0});
    affine<float> translation(2,params3);

    translation.print("translation");
    translation.print_data();

    x1 = translation.apply(x0);
    x1.print_data("grid x1, translation [10, -10]");

    return 0;
};