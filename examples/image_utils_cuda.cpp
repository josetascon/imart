/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-28 10:02:24
*/

// File to test utilities such as: pad, unpad, normalize

// std libs
#include <iostream>
#include <memory>
#include <complex>

// local libs
#include "../src/image.h"

using namespace imart;
using type = float;

int main()
{
    // ============================================
    //      Testing image utility functions
    // ============================================
    // test padding
    image_cuda<type> img10(8,5);
    img10.ones();
    img10.random();
    image_cuda<type> img10pad = pad(img10, std::vector<int>{1,3}, std::vector<int>{4,2});
    
    std::cout << "===================== ";
    std::cout << "Test function pad";
    std::cout << " =====================";
    std::cout << std::endl;
    std::cout << "pad width (1 column left, 4 column right)\n";
    std::cout << "pad height (3 rows up, 2 rows down)\n";
    img10.print("img info");
    img10.print_data();
    img10pad.print_data("img to pad");
    // std::cout << std::endl;

    std::cout << "===================== ";
    std::cout << "Test function unpad";
    std::cout << " =====================";
    std::cout << std::endl;
    
    image_cuda<type> img20 = unpad(img10pad, std::vector<int>{1,3}, std::vector<int>{4,2});
    // img20.print("img unpad");
    img20.print_data("unpad to get img");

    // ============================================
    //      Testing image 3d
    // ============================================
    // test padding
    image_cuda<type> img13(4,3,2);
    img13.ones();
    img13.random();
    image_cuda<type> img13pad = pad(img13, std::vector<int>{1,3,2}, std::vector<int>{4,2,1});
    
    std::cout << "===================== ";
    std::cout << "Test function pad 3d";
    std::cout << " =====================";
    std::cout << std::endl;
    std::cout << "pad width (1 column left, 4 column right)\n";
    std::cout << "pad height (3 rows up, 2 rows down)\n";
    std::cout << "pad depth (2 depth up, 1 depth down)\n";
    img13.print("img info");
    img13.print_data();
    img13pad.print_data("img to pad");
    // std::cout << std::endl;

    std::cout << "===================== ";
    std::cout << "Test function unpad 3d";
    std::cout << " =====================";
    std::cout << std::endl;
    
    image_cuda<type> img23 = unpad(img13pad, std::vector<int>{1,3,2}, std::vector<int>{4,2,1});
    // img20.print("img unpad");
    img23.print_data("unpad to get img");


    std::cout << "===================== ";
    std::cout << "Test function cast";
    std::cout << " =====================";
    std::cout << std::endl;
    image_cuda<float> img30(4,3);
    image_cuda<int> img31;
    image_cuda<double> img32;
    img30.random(0.0, 100.0);

    // cast
    cast(img30, img31);
    cast(img30, img32);
    
    img30.print_data("float");
    img31.print_data("cast to int");
    img32.print_data("cast to double");
    
    return 0;
};