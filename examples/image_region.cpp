/*
* @Author: Jose Tascon
* @Date:   2019-11-07 10:13:08
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-07-15 16:23:46
*/


// std libs
#include <iostream>

// local libs
#include "../src/image.h"

using namespace imart;

int main()
{
    using type = float;

    // 2d
    std::cout<< "========================================" << std::endl;
    std::cout<< "          Two dimensional Region" << std::endl;
    std::cout<< "========================================" << std::endl;
    auto img10 = image_cpu<type>::new_pointer(8,5);
    img10->random();
    img10->print();
    img10->print_data("Random values");
    
    auto img11 = img10->region(std::vector<int>{2,1}, std::vector<int>{4,2});
    img11->print_data("Extract Region");

    // 3d
    // std::cout<< "========================================" << std::endl;
    // std::cout<< "          Three dimensional Region" << std::endl;
    // std::cout<< "========================================" << std::endl;
    

    return 0;
}