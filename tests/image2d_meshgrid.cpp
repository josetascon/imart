/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-01-24 16:32:16
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image_2d.h"
#include "../src/image_3d.h"
#include "../src/grid.h"

int main()
{
    // ============================================
    //      Testing image_2d meshgrid
    // ============================================
    // Create small imame
    image<double> image0(10,5);

    grid<double> x0;
    std::cout << "a\n" ;
    x0.print();
    
    x0.meshgrid(image0);
    std::cout << "b\n" ;

    x0.print();
    x0.print_data();

    // ============================================
    //      Testing image_3d meshgrid
    // ============================================
    image<double> image1(3, 5, 2);
    grid<double> x1;
    std::cout << "a\n" ;
    x1.print();
    
    x1.meshgrid(image1);
    std::cout << "b\n" ;

    x1.print();
    x1.print_data();

    // std::cout << "Image(10,5) mesh grid:" << std::endl;
    // std::cout << "x:" << std::endl;
    // (*image0.get_grid())[0].print_data();
    // std::cout << "y:" << std::endl;
    // (*image0.get_grid())[1].print_data();

    // // Create small imame
    // image_2d<float> image1(1000,200);

    // image1.meshgrid();
    // // image1.grid[0].print_data();

    // // Create small imame
    // image_2d<float> image2(1000,200);

    // image2.meshgrid();

    return 0;

};