/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2019-12-04 11:34:42
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image_2d.h"
#include "../src/grid.h"
#include "../src/affine_2d.h"
#include "../src/interpolate.h"

// using type = int;
using type = float;
// using type = double;

image_2d<type> ssd(image_2d<type> & image0, image_2d<type> & image1, affine_2d<type> & transform)
{
    // image0.print_data();
    // image1.print_data();

    grid<type> x0(image0);
    grid<type> x1(image1);
    interpolate<type> image1_p(image1, x1);

    // image_2d<type> ssd_(image0.get_width(),image0.get_height());
    // ssd_ = (image1_p*(transform * x0));
    image_2d<type> ssd_ = image0 - image1_p*(transform * x0);

    return ssd_;
};


int main()
{
    // ============================================
    //          Testing interpolate
    // ============================================

    // Create small imame
    // image_2d<type> image0(4,3);
    image_2d<type> image0(10,7);
    // image0.random();
    image0.ones();
    grid<type> x0(image0);

    image_2d<type> params(6,1);
    *params.get_data() = {1.0, 0.0, 0.0, 1.0, 3.5, -1.0};
    // std::shared_ptr<std::vector<type>> buffer = std::make_shared<std::vector<type>>(6);
    // *buffer = {1.0, 0.0, 0.0, 1.0, 1.5, 1.0};
    // *buffer = {1, 0, 0, 1, 4, 0};
    // image_2d<type> params(buffer, 6,1);
    affine_2d<type> translation(params);

    grid<type> t_x0 = translation.transform(x0);
    // x0.print_data();
    t_x0.print_data();
    translation.print();

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


    // Testing the inverse transfor applied to ssd function
    std::cout << "Test compressed notation\n";
    std::shared_ptr<std::vector<type>> buffer_inv = std::make_shared<std::vector<type>>(6);
    *buffer_inv = {1.0, 0.0, 0.0, 1.0, -1.5, -1.0};
    image_2d<type> params_inv(buffer_inv, 6,1);
    affine_2d<type> inv_translation(params_inv);

    image_2d<type> result = ssd(image0, image0t, inv_translation);
    result.print_data("SSD now as image:");
    // image0.print_data();
    // image0t.print_data();

    // image0t.print_data();

    // ssd metric

    // type ssd(image0 , image1, transform)
    // {
    //     grid<type> x0(image0);
    //     grid<type> x1(image1);
    //     interpolate<type> image1_p(image1, x1)

    //     image_2d<type> result = (image0 - image1_p*transform*x0)^2.0;
    //     return result.sum()
    // };

    


    return 0;

};