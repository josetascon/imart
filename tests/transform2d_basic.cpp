/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2019-11-19 00:13:42
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image_2d.h"
#include "../src/affine_2d.h"

int main()
{
    // ============================================
    //          Testing affine_2d 
    // ============================================
    // Create small imame
    std::shared_ptr<float[]> buffer(new float[6] {1.1, 2.1, 3.1, 4.1, 5.1, 6.1});
    image_2d<float> params(buffer, 6,1);
    affine_2d<float> affine1(params);

    std::vector<float> point1({1.0,2.0});
    std::vector<float> point2;

    // affine1.print();
    std::cout << affine1;

    point2 = affine1.transform(point1);
    std::cout << point1[0] << " " << point1[1] << std::endl;
    std::cout << point2[0] << " " << point2[1] << std::endl;


    // image_2d<float> image1(4,3);
    // image1.random();

    // affine1.transform(image1);


    // affine1





};