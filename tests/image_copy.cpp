/*
* @Author: Jose Tascon
* @Date:   2020-01-27 15:51:36
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-02-06 15:35:31
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/transform_base.h"

using type = float;

int main()
{
    std::cout << "===================== ";
    std::cout << "Test image, equal and duplicate functions";
    std::cout << " =====================";
    std::cout << std::endl;

    image<type> image1(6,1);
    image<type> image2(6,1);
    image<type> image3(6,1);
    *image1.get_data() = {1.0, 0.0, 0.0, 1.0, 3.5, -1.0};

    image2 = image1;
    image3.duplicate(image1);

    std::cout << "Summary commands:\n";
    std::cout << "image1 = {1.0, 0.0, 0.0, 1.0, 3.5, -1.0}\n";
    std::cout << "image2 = image1\n";
    std::cout << "image3.duplicate(image1)\n";
    std::cout << std::endl;

    image1.print_data("Image1:");
    image2.print_data("Image2:");
    image3.print_data("Image3:");

    std::cout << "Image 1 ptr count: " << image1.get_ptr_count() << "\n";
    std::cout << "Image 1 ptr value: " << image1.ptr() << "\n";
    std::cout << "Image 2 ptr value: " << image2.ptr() << "\n";
    std::cout << "Image 3 ptr value: " << image3.ptr() << "\n";
    if (image1.ptr() == image2.ptr() and image1.ptr() == image3.ptr()) { std::cout << "Same pointers!\n"; };
    std::cout << std::endl;

    std::cout << "===================== ";
    std::cout << "Test image, copy functions";
    std::cout << " =====================";
    std::cout << std::endl;

    image<type> image11(6,1);
    image<type> image12(6,1);
    *image11.get_data() = {9.0, 0.0, 0.0, 9.0, 3.5, -1.0};
    image12.copy(image11);

    std::cout << "Summary commands:\n";
    std::cout << "image11 = {9.0, 0.0, 0.0, 9.0, 3.5, -1.0}\n";
    std::cout << "image12.copy(image11)\n";
    std::cout << std::endl;

    image11.print_data("Image11:");
    image12.print_data("Image12:");

    std::cout << "Image 11 ptr count: " << image11.get_ptr_count() << "\n";
    std::cout << "Image 11 ptr value: " << image11.ptr() << "\n";
    std::cout << "Image 12 ptr value: " << image12.ptr() << "\n";
    if (image11.ptr() != image12.ptr()) { std::cout << "Different pointers!\n"; };


    std::cout << "===================== ";
    std::cout << "Test transform, equal and duplicate functions";
    std::cout << " =====================";
    std::cout << std::endl;

    
    image<type>::pointer imaget1(new image<type>{1.0, 0.0, 0.0, 1.0, 0.5, -0.5});
    
    transform_base<type> transform1(2,imaget1);
    transform_base<type> transform2;
    transform_base<type> transform3;

    transform2 = transform1;
    transform3.duplicate(transform1);

    std::cout << "Summary commands:\n";
    std::cout << "transform1 with params = {1.0, 0.0, 0.0, 1.0, 3.5, -1.0}\n";
    std::cout << "transform2 = transform1\n";
    std::cout << "transform3.duplicate(transform1)\n";
    std::cout << std::endl;

    transform1.print_data("transform1:");
    transform2.print_data("transform2:");
    transform3.print_data("transform3:");

    // image<type>::pointer imaget2 = transform1.get_parameters();
    // std::cout << "Ptr: " << imaget1->ptr() << std::endl;
    // std::cout << "Ptr: " << imaget2->ptr() << std::endl;
    // std::cout << "Ptr count: " << imaget1->get_ptr_count() << std::endl;
    // std::cout << "Ptr count: " << imaget2->get_ptr_count() << std::endl;

    image<type>::pointer param_t1 = transform1.get_parameters();
    image<type>::pointer param_t2 = transform2.get_parameters();
    image<type>::pointer param_t3 = transform3.get_parameters();

    std::cout << "Transform 1 parameters (image) ptr count: " << param_t1.use_count() << "\n";    
    std::cout << "Transform 1 parameters (image) ptr value: " << param_t1.get() << "\n";
    std::cout << "Transform 2 parameters (image) ptr value: " << param_t2.get() << "\n";
    std::cout << "Transform 3 parameters (image) ptr value: " << param_t3.get() << "\n";

    // The transform is duplicating the image pointer. Not the inner image data pointer
    std::cout << "Transform 1 parameters (image data) ptr count: " << imaget1->get_ptr_count() << "\n";
    std::cout << "Transform 1 parameters (image data) ptr value: " << param_t1->ptr() << "\n";
    std::cout << "Transform 2 parameters (image data) ptr value: " << param_t2->ptr() << "\n";
    std::cout << "Transform 3 parameters (image data) ptr value: " << param_t3->ptr() << "\n";
    if (param_t1->ptr() == param_t2->ptr() and param_t1->ptr() == param_t3->ptr()) { std::cout << "Same pointers!\n"; };
    std::cout << std::endl;

    std::cout << "===================== ";
    std::cout << "Test transform, copy functions";
    std::cout << " =====================";
    std::cout << std::endl;

    image<type>::pointer imaget11(new image<type>{9.0, 0.0, 0.0, 9.0, 3.5, -1.0});
    transform_base<type> transform11(2,imaget11);
    transform_base<type> transform12(2);
    
    transform12.copy(transform11);

    std::cout << "Summary commands:\n";
    std::cout << "transform11 with params = {9.0, 0.0, 0.0, 9.0, 3.5, -1.0}\n";
    std::cout << "transform12.copy(transform11)\n";
    std::cout << std::endl;
    
    transform11.print_data("transform11:");
    transform12.print_data("transform12:");

    std::cout << "transform11 ptr count: " << (transform11.get_parameters())->get_ptr_count() << "\n";
    std::cout << "transform11 ptr value: " << (transform11.get_parameters())->ptr() << "\n";
    std::cout << "transform12 ptr value: " << (transform12.get_parameters())->ptr() << "\n";
    if (image11.ptr() != image12.ptr()) { std::cout << "Different pointers!\n"; };


    return 0;

};