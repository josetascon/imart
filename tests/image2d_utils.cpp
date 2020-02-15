/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-02-12 00:46:10
*/

// File to test utilities such as: pad, unpad, normalize

// std libs
#include <iostream>
#include <memory>
#include <complex>

// local libs
#include "../src/image.h"

using type = float;

int main()
{
    // ============================================
    //      Testing image utility functions
    // ============================================
    // // test padding
    image<type> img10(8,5);
    img10.ones();
    img10.random();
    image<type> img10pad = pad(img10, std::vector<int>{1,3}, std::vector<int>{4,2});
    
    std::cout << "===================== ";
    std::cout << "Test function pad";
    std::cout << " =====================";
    std::cout << std::endl;
    std::cout << "pad width (1 column left, 4 column right)\n";
    std::cout << "pad height (3 rows up, 2 rows down)\n";
    img10.print("im info");
    img10.print_data();
    img10pad.print_data("im to pad");
    // std::cout << std::endl;

    // std::cout << "===================== ";
    // std::cout << "Test function unpad";
    // std::cout << " =====================";
    // std::cout << std::endl;
    
    image<type> img20 = unpad(img10pad, std::vector<int>{1,3}, std::vector<int>{4,2});
    // img20.print("im unpad");
    img20.print_data("unpad to get im");


    std::cout << "===================== ";
    std::cout << "Test function cast";
    std::cout << " =====================";
    std::cout << std::endl;
    image<int> img30;
    image<float> img31(4,3);
    image<double> img32;
    int v1;
    double v2;

    img31.random(0.0, 100.0);

    img30 = cast(img31, v1);
    img32 = cast(img31, v2);
    
    img31.print_data("float");
    img30.print_data("cast to int");
    img32.print_data("cast to double");



    return 0;
};