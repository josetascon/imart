/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2019-11-11 15:08:23
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../inc/image_base.hpp"

int main()
{
    // ============================================
    //      Testing ImageBase2D basic features
    // ============================================
    // Create empty imame
    image_base_2d<int> image0;     

    std::cout << "===================== ";
    std::cout << "Test class image_base_2d, basic features";
    std::cout << " =====================";
    std::cout << std::endl;

    std::cout << "Creating an empty image";
    std::cout << std::endl;    

    std::cout << "Image width: ";
    std::cout << image0.get_width();
    std::cout << std::endl;

    std::cout << "Image height: ";
    std::cout << image0.get_height();
    std::cout << std::endl;

    std::cout << "Image elements: ";
    std::cout << image0.get_total_elements();
    std::cout << std::endl;

    std::cout << "Image data pointer: ";
    std::cout << image0.get_data();
    std::cout << std::endl;
    std::cout << "internal image ptr count: ";
    std::cout << image0.get_ptr_count();
    std::cout << std::endl;
    std::cout << "external image ptr count: ";
    std::cout << image0.get_data().use_count();
    std::cout << std::endl;


    // Create medium size image 
    image_base_2d<float> image1(4, 3);

    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image_base_2d, small size image";
    std::cout << " =====================";
    std::cout << std::endl;

    std::cout << "Creating an image with (4, 3)";
    std::cout << std::endl;    

    std::cout << "Image width: ";
    std::cout << image1.get_width();
    std::cout << std::endl;

    std::cout << "Image height: ";
    std::cout << image1.get_height();
    std::cout << std::endl;

    std::cout << "Image elements: ";
    std::cout << image1.get_total_elements();
    std::cout << std::endl;

    std::cout << "Image data pointer: ";
    std::cout << image1.get_data();
    std::cout << std::endl;
    image1.print_ptr_count();   // testing print function, 
                                // equivalent to the commented lines
    // std::cout << "internal image ptr count: ";
    // std::cout << image1.get_ptr_count();
    // std::cout << std::endl;
    std::cout << "external image ptr count: ";
    std::cout << image1.get_data().use_count();
    std::cout << std::endl;

    image1.print_data();


    // Create medium size image 
    image_base_2d<double> image2(128, 64);

    // std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image_base_2d, medium size image";
    std::cout << " =====================";
    std::cout << std::endl;

    std::cout << "Creating an image with (128, 64)";
    std::cout << std::endl;    

    std::cout << "Image width: ";
    std::cout << image2.get_width();
    std::cout << std::endl;

    std::cout << "Image height: ";
    std::cout << image2.get_height();
    std::cout << std::endl;

    std::cout << "Image elements: ";
    std::cout << image2.get_total_elements();
    std::cout << std::endl;

    std::cout << "Image data pointer: ";
    std::cout << image2.get_data();
    std::cout << std::endl;
    std::cout << "internal image ptr count: ";
    std::cout << image2.get_ptr_count();
    std::cout << std::endl;
    std::cout << "external image ptr count: ";
    std::cout << image2.get_data().use_count();
    std::cout << std::endl;
    

    // Create different image objects
    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image_base_2d, create and equal";
    std::cout << " =====================";
    std::cout << std::endl;
    
    std::shared_ptr<float[]> buffer(new float[6] {1.1, 2.1, 3.1, 4.1, 5.1, 6.1});

    image_base_2d<float> image3(5,3);
    image_base_2d<float> image4(buffer, 3,2);
    

    std::cout << "image3 ptr count: " << image3.get_ptr_count() << std::endl;
    std::cout << "image4 ptr count: " << image4.get_ptr_count() << std::endl;

    image3.print("image3");
    image4.print("image4");

    image3 = image4;
    std::cout << std::endl << "Making image3 = image4. " << std::endl;

    image3.print("image3");
    image4.print("image4");
    std::cout << std::endl;

    image3.print_data();

    std::cout << "image3 ptr count: " << image3.get_ptr_count() << std::endl;
    std::cout << "image4 ptr count: " << image4.get_ptr_count() << std::endl;
    

    return 0;
};