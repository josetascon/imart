/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-18 06:34:35
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/grid.h"

using namespace imart;

int main()
{
    // ============================================
    //      Testing image_2d meshgrid
    // ============================================
    // Create small imame
    image_gpu<double> image0(10,5);
    
    grid<double,vector_gpu<double>> x0(image0);
    std::cout << "a\n" ;
    x0.print();
    
    x0.meshgrid();
    std::cout << "b\n" ;
    x0.print_data();

    image_gpu<double> image1(3, 5, 2);
    grid_gpu<double> x1; // default grid in 2d
    std::cout << "a\n" ;
    x1.print();
    
    x1.meshgrid(image1);
    std::cout << "b\n" ;

    x1.print();
    x1.print_data();

    return 0;

};