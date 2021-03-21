/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-10 16:48:56
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
    // Create empty image
    image_opencl<int> image0;

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
    image_opencl<float> image1(4, 3);

    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, small size image";
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

    // std::cout << "Image data pointer: ";
    // std::cout << image1.get_data();
    // std::cout << std::endl;
    // image1.print_ptr_count();   // testing print function, 
    //                             // equivalent to the commented lines
    // // std::cout << "internal image ptr count: ";
    // // std::cout << image1.get_ptr_count();
    // // std::cout << std::endl;
    // std::cout << "external image ptr count: ";
    // std::cout << image1.get_data().use_count();
    // std::cout << std::endl;

    // image1.print_data();

    /*
    // Create medium size image 
    image<double> image2(128, 64);

    // std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, medium size image";
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
    */

    // Create different image objects
    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, create and equal";
    std::cout << " =====================";
    std::cout << std::endl;
    
    std::shared_ptr<std::vector<float>> buffer = std::make_shared<std::vector<float>>(6);
    *buffer = {1.1, 2.1, 3.1, 4.1, 5.1, 6.1};

    // image_opencl<float> image3(5,3);
    image_opencl<float> image4(buffer, 3,2);
    // image4 = image4 + image4;

    // std::cout << "image3 ptr count: " << image3.get_ptr_count() << std::endl;
    // std::cout << "image4 ptr count: " << image4.get_ptr_count() << std::endl;

    // image3.print("image3");
    image4.print("image4");
    image4.print_data();

    /*
    image3 = image4;
    std::cout << std::endl << "Making image3 = image4. " << std::endl;

    image3.print("image3");
    image4.print("image4");
    std::cout << std::endl;

    image3.print_data();

    std::cout << "image3 ptr count: " << image3.get_ptr_count() << std::endl;
    std::cout << "image4 ptr count: " << image4.get_ptr_count() << std::endl;

    std::cout << std::endl << "\nMaking new image-pointer (2,2)\n";  
    image<float>::pointer image_ptr( new image<float>(2,2) );
    // image<float>::pointer image_ptr2 = image<float>::new_pointer(); // ** to implement in future
    
    image_ptr->print("image-pointer");
    image_ptr->print_data();

    std::cout << "image+pointer ptr count: " << image_ptr->get_ptr_count() << std::endl;
    
    *image_ptr = image3;
    std::cout << std::endl << "Making image-pointer = image3. " << std::endl;

    image_ptr->print("image-pointer");
    std::cout << "image3 ptr count: " << image3.get_ptr_count() << std::endl;
    std::cout << "image-pointer ptr count: " << image_ptr->get_ptr_count() << std::endl;

    std::cout << std::endl << "Testing constructor with image(input.get_size()). " << std::endl;
    image<double> image_size(image3.get_size());
    image_size.print("image construted with image3.get_size()");

    // Create different image objects
    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, initilization";
    std::cout << " =====================";
    std::cout << std::endl;

    image<float> image_list{1.2, 2.1, 1.3, 0.0, 3.1, 1.5, 5.1};
    std::cout << "image initialized with list {1.2, 2.1, 1.3, 0.0, 3.1, 1.5, 5.1}" << std::endl;
    image_list.print("image list");
    image_list.print_data();

    image<float> image5(6,6);
    // image<unsigned int> image6(4,7);
    
    image5.print("image5");
    image5.zeros();
    image5.print_data("zeros:");
    image5.ones();
    image5.print_data("ones:");
    image5.random();
    image5.print_data("random:");


    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, access elements";
    std::cout << " =====================";
    std::cout << std::endl;
    std::cout << "image5(0) = " << image5(0) << std::endl;
    std::cout << "image5(11) = " << image5(11) << std::endl;
    std::cout << "image5(35) = " << image5(35) << std::endl;
    std::cout << "image5(0,0) = " << image5(0,0) << std::endl;
    std::cout << "image5(2,3) = " << image5(2,3) << std::endl;
    std::cout << "image5(5,4) = " << image5(5,4) << std::endl;
    std::cout << "image5(5,5) = " << image5(5,5) << std::endl;
    */
    return 0;
};