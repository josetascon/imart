/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-02-11 20:09:33
*/

// std libs
#include <iostream>
#include <memory>
#include <complex>

// local libs
#include "../src/image.h"

int main()
{
    // ============================================
    //      Testing image basic features
    // ============================================
    // Create empty image
    image<std::complex<float>> image1(6,6);
    image<std::complex<float>> image2(6,6);
    image<std::complex<float>> image3(6,6);
    
    std::cout << "===================== ";
    std::cout << "Test class image, initialize";
    std::cout << " =====================";
    image1.print("image1");
    image1.zeros();
    image1.print_data("zeros:");
    image1.ones();
    image1.print_data("ones:");
    image1.random();
    image1.print_data("random:");



    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, access elements";
    std::cout << " =====================";
    std::cout << std::endl;
    std::cout << "image1(0) = " << image1(0) << std::endl;
    std::cout << "image1(11) = " << image1(11) << std::endl;
    std::cout << "image1(35) = " << image1(35) << std::endl;
    std::cout << "image1(0,0) = " << image1(0,0) << std::endl;
    std::cout << "image1(2,3) = " << image1(2,3) << std::endl;
    std::cout << "image1(5,4) = " << image1(5,4) << std::endl;
    std::cout << "image1(5,5) = " << image1(5,5) << std::endl;
    std::cout << std::endl;



    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, complex operations";
    std::cout << " =====================";
    std::cout << std::endl;

    image1.random();
    image1.print_data("image 1 with random values:");
    image2.fill(std::complex<float>(1.0,-0.3));
    image2.print_data("image2 with same values:");

    std::cout << "image3 = image1 + image2" << std::endl;

    image3 = image1+image2;
    image3.print_data();

    return 0;
};