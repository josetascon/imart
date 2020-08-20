/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-19 16:23:46
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"
#include "../src/inearest.h"

using namespace imart;

using type = float;
// using type = double;

int main()
{
    // ============================================
    //   Testing interpolator after transformation
    // ============================================

    // Create small image and grid
    image_cpu<type>::pointer image0 = image_cpu<type>::new_pointer(7,5);
    grid_cpu<type>::pointer x0 = grid_cpu<type>::new_pointer(*image0);

    // Initialize image
    image0->random();
    // image0->ones();
    image0->print_data("image 0:");
    x0->print_data("grid x0:");
    
    // Create affine transform
    image_cpu<type>::pointer params( new image_cpu<type>{1.0, 0.0, 0.0, 1.0, 2.5, -1.0} );
    affine<type> translation( 2, params );

    // Apply transform
    grid_cpu<type>::pointer x1 = grid_cpu<type>::new_pointer();
    x1 = translation.apply(x0);
    translation.print_data("Transform parameters:");
    x1->print_data("x1 = transform(x0)");

    // Create interpolator
    // inearest<type,vector_cpu<type>> image0_intp(image0, x0);
    inearest<type,vector_cpu<type>> image0_intp(image0); // new version
    image_cpu<type>::pointer image1;
    image1 = image0_intp.apply(x1);
    image1->print_data("image 1:");

    // ============================================
    //   Testing interpolator with manual grid
    // ============================================
    // Create a grid from scratch
    // grid_cpu<type>::pointer x2 = grid_cpu<type>::new_pointer(2);
    // x2->set_size(std::vector<int>{7, 5});
    grid_cpu<type>::pointer x2 = grid_cpu<type>::new_pointer(std::vector<int>{7, 5});
    x2->set_spacing(std::vector<double>{2.0, 1.5});
    x2->set_origin(std::vector<double>{-1.0, 2.5});
    x2->meshgrid();
    x2->print_data("grid x2:");

    // Interpolate
    image_cpu<type>::pointer image2;
    image2 = image0_intp.apply(x2);
    image2->print_data("image 2:");


    return 0;

};