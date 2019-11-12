/*
* @Author: Jose Tascon
* @Date:   2019-11-07 10:13:08
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2019-11-12 09:56:23
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../inc/image_base.hpp"

int main()
{
    // ============================================
    //      Testing ImageBase2D basic operations
    // ============================================

    std::shared_ptr<float[]> buffer1(new float[6] {1.1, 2.1, 3.1, 4.1, 5.1, 6.1});
    std::shared_ptr<float[]> buffer2(new float[6] {2.1, 1.1, 0.1, -0.9, -1.9, -2.9});
    

    image_base_2d<float> image1(buffer1, 3,2);
    image_base_2d<float> image2(buffer2, 3,2);
    image_base_2d<float> image3;

    std::cout << "Adding 2 images:" << std::endl;
    std::cout << "image3 = image1 + image2" << std::endl;
    image3 = image1 + image2;

    image3.print("image3 Info");
    
    image1.print_data("image1:");
    image2.print_data("image2:");
    image3.print_data("image3:");

    std::cout << "image1 ptr count: " << image1.get_ptr_count() << std::endl;
    std::cout << "image2 ptr count: " << image2.get_ptr_count() << std::endl;
    std::cout << "image3 ptr count: " << image3.get_ptr_count() << std::endl;



    return 0;
}