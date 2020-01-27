/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-01-27 09:26:35
*/

// std libs
#include <iostream>
#include <memory>
#include <vector>

// local libs
#include "../src/image.h"
#include "../src/affine.h"

int main()
{
    // ============================================
    //      Testing affine with points
    // ============================================
    // Create identity transform
    std::shared_ptr<std::vector<float>> buffer1 = std::make_shared<std::vector<float>>(6);
    *buffer1 = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    image<float> params1(buffer1, 6,1);
    affine<float> affine1(2,params1);

    std::vector<float> point1({1.0,2.0});
    std::vector<float> point2;

    // affine1.print();
    std::cout << affine1;

    point2 = affine1.transform(point1);
    std::cout << "Input point:" << std::endl;
    std::cout << point1[0] << " " << point1[1] << std::endl;
    std::cout << "Transformed point (identity):" << std::endl;
    std::cout << point2[0] << " " << point2[1] << std::endl;

    // Create other transform
    std::shared_ptr<std::vector<float>> buffer2 = std::make_shared<std::vector<float>>(6);
    *buffer2 = {1.1, 0.5, -0.5, 0.9, 2.1, -1.1};
    image<float> params2(buffer2, 6,1);
    affine<float> affine2(2,params2);

    std::vector<float> point3({1.0,2.0});
    std::vector<float> point4;

    std::cout << affine2;

    point4 = affine2.transform(point3);
    std::cout << "Input point:" << std::endl;
    std::cout << point3[0] << " " << point3[1] << std::endl;
    std::cout << "Transformed point:" << std::endl;
    std::cout << point4[0] << " " << point4[1] << std::endl;
    std::cout << std::endl;
    
    // ============================================
    //      Testing affine with grids
    // ============================================
    image<float> image0(5,3);
    grid<float> x0(image0);
    grid<float> x1;

    // x0.print();
    x1 = affine1.transform(x0);

    x0.print_data("grid x0");
    x1.print_data("grid x1, affine 1 (identity)");

    x1 = affine2.transform(x0);
    x1.print_data("grid x1, affine 2");

    std::shared_ptr<std::vector<float>> buffer3 = std::make_shared<std::vector<float>>(6);
    *buffer3 = {1.0, 0.0, 0.0, 1.0, 10.0, -10.0};
    image<float> params3(buffer3, 6,1);
    affine<float> translation(2,params3);

    x1 = translation.transform(x0);
    x1.print_data("grid x1, translation [10, -10]");    

    return 0;
};