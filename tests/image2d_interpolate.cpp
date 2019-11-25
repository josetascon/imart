/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2019-11-22 16:43:00
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image_2d.h"
#include "../src/grid.h"
#include "../src/affine_2d.h"
#include "../src/interpolate.h"


int main()
{
    // ============================================
    //          Testing interpolate
    // ============================================

    // using type = int;
    using type = float;
    // using type = double;


    // Create small imame
    image_2d<type> image0(10,7);
    image0.random();
    grid<type> x0(image0);

    std::shared_ptr<std::vector<type>> buffer = std::make_shared<std::vector<type>>(6);
    *buffer = {1.0, 0.0, 0.0, 1.0, 4.0, 0.0};
    // *buffer = {1, 0, 0, 1, 4, 0};
    image_2d<type> params(buffer, 6,1);
    affine_2d<type> translation(params);

    grid<type> t_x0 = translation.transform(x0);
    // x0.print_data();
    t_x0.print_data();

    std::cout << "Create interpolation object\n";
    interpolate<type> image1_itp(image0, x0);
    std::cout << "Linear interpolation function\n";
    image_2d<type> image0t = image1_itp.linear(t_x0); 
    // TODO: Error of free() when using double. The error is in this method
    // *** TEST neighbors4??? change array to vector
    

    std::cout << "Linear interpolation finished\n";

    image0.print();
    image0.print_data();
    
    image0t.print();
    image0t.print_data();


    // ssd metric

    // float ssd(image0 , image1, transform)
    // {
    //     grid<float> x0(image0);
    //     grid<float> x1(image1);
    //     interpolate<float> image1_p(image1, x1)

    //     image_2d<float> result = (image0 - image1_p*transform*x0)^2.0;
    //     return result.sum()
    // };

    


    return 0;

};