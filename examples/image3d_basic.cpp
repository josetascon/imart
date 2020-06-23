/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-02-15 16:08:48
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"

int main()
{
    // ============================================
    //      Testing image basic features
    // ============================================
    // Create empty imame
    image<int> image0(3);

    std::cout << "===================== ";
    std::cout << "Test class image, basic features";
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

    std::cout << "Image length: ";
    std::cout << image0.get_length();
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

    std::cout << image0;


    // Create medium size image 
    image<float> image1(4, 3, 2);

    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, small size image";
    std::cout << " =====================";
    std::cout << std::endl;

    std::cout << "Creating an image with (4, 3, 2)";
    std::cout << std::endl;    

    std::cout << "Image width: ";
    std::cout << image1.get_width();
    std::cout << std::endl;

    std::cout << "Image height: ";
    std::cout << image1.get_height();
    std::cout << std::endl;

    std::cout << "Image length: ";
    std::cout << image1.get_length();
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
    image<double> image2(128, 64, 32);

    // std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, medium size image";
    std::cout << " =====================";
    std::cout << std::endl;

    std::cout << "Creating an image with (16, 64, 32)";
    std::cout << std::endl;    

    std::cout << "Image width: ";
    std::cout << image2.get_width();
    std::cout << std::endl;

    std::cout << "Image height: ";
    std::cout << image2.get_height();
    std::cout << std::endl;

    std::cout << "Image length: ";
    std::cout << image2.get_length();
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
    std::cout << "Test class image, create and equal";
    std::cout << " =====================";
    std::cout << std::endl;
    
    
    std::shared_ptr<std::vector<float>> buffer = std::make_shared<std::vector<float>>(12);
    *buffer = {1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 6.1, 5.1, 4.1, 3.1, 2.1, 1.1};

    image<float> image3(5,3,4);
    image<float> image4(buffer, 3,2,2);
    // image4 = image4 + image4;

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

    std::cout << std::endl << "Testing constructor with image(input.get_size()). " << std::endl;
    image<double> image_size(image3.get_size());
    image_size.print("image construted with image3.get_size()");

    // Create different image objects
    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, initilization";
    std::cout << " =====================";
    std::cout << std::endl;

    image<float> image5(4,4,4);
    // image<unsigned int> image6(4,7);
    
    image5.print("image5");
    image5.zeros();
    image5.print_data("zeros:");
    image5.ones();
    image5.print_data("ones:");
    image5.random();
    // image5(0) = 200.1;
    // float value = image5(0);
    // value = 300.1;
    image5.print_data("random:");


    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, access elements";
    std::cout << " =====================";
    std::cout << std::endl;
    std::cout << "image5(0) = " << image5(0) << std::endl;
    std::cout << "image5(31) = " << image5(31) << std::endl;
    std::cout << "image5(63) = " << image5(63) << std::endl;
    std::cout << "image5(0,0,0) = " << image5(0,0,0) << std::endl;
    std::cout << "image5(0,2,1) = " << image5(0,2,1) << std::endl;
    std::cout << "image5(3,1,2) = " << image5(3,1,2) << std::endl;
    std::cout << "image5(3,3,3) = " << image5(3,3,3) << std::endl;
    
    return 0;
};