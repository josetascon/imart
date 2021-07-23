/*
* @Author: Jose Tascon
* @Date:   2019-11-07 10:13:08
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-07-17 06:10:51
*/


// std libs
#include <iostream>

// local libs
#include "../src/image.h"
#include "../src/image_utils.h"

using namespace imart;

int main()
{
    using type = unsigned char;

    // 2d
    std::cout<< "========================================" << std::endl;
    std::cout<< "       Two dimensional Bounding Box" << std::endl;
    std::cout<< "========================================" << std::endl;
    auto img10 = image_cpu<type>::new_pointer();
    img10->read("examples/images/t0.png");
    std::vector<std::vector<int>> bbox = bounding_box(img10);

    std::cout << "x: " << bbox[0][0] << std::endl;
    std::cout << "y: " << bbox[0][1] << std::endl;
    std::cout << "w: " << bbox[1][0] << std::endl;
    std::cout << "h: " << bbox[1][1] << std::endl;


    // 3d
    // std::cout<< "========================================" << std::endl;
    // std::cout<< "          Three dimensional Region" << std::endl;
    // std::cout<< "========================================" << std::endl;
    

    return 0;
}