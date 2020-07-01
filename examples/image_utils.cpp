/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-26 16:09:05
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
    //      Testing image 2d utility functions
    // ============================================
    // test padding
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
    img10.print("img info");
    img10.print_data();
    img10pad.print_data("img to pad");
    // std::cout << std::endl;

    std::cout << "===================== ";
    std::cout << "Test function unpad";
    std::cout << " =====================";
    std::cout << std::endl;
    
    image<type> img20 = unpad(img10pad, std::vector<int>{1,3}, std::vector<int>{4,2});
    // img20.print("img unpad");
    img20.print_data("unpad to get img");

    // ============================================
    //      Testing image 3d
    // ============================================
    // test padding
    image<type> img13(4,3,2);
    img13.ones();
    img13.random();
    image<type> img13pad = pad(img13, std::vector<int>{1,3,2}, std::vector<int>{4,2,1});
    
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
    
    image<type> img23 = unpad(img13pad, std::vector<int>{1,3,2}, std::vector<int>{4,2,1});
    // img20.print("img unpad");
    img23.print_data("unpad to get img");


    std::cout << "===================== ";
    std::cout << "Test function cast";
    std::cout << " =====================";
    std::cout << std::endl;
    image<int> img30;
    image<float> img31(4,3);
    image<double> img32;
    // int v1;
    // double v2;

    img31.random(0.0, 100.0);

    // img30 = cast<int,vector_cpu<float>,float>(img31);
    img30 = cast<int,vector_cpu<int>,float,vector_cpu<float>>(img31);
    img32 = cast<double,vector_cpu<double>,float,vector_cpu<float>>(img31);
    
    img31.print_data("float");
    img30.print_data("cast to int");
    img32.print_data("cast to double");



    return 0;
};