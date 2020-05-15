/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-03-23 14:29:50
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"
#include "../src/interpolate.h"

// using type = int;
using type = float;
// using type = double;

// image<type> ssd(image<type> & image0, image<type> & image1, affine<type> & transform)
// {
//     // image0.print_data();
//     // image1.print_data();

//     grid<type> x0(image0);
//     grid<type> x1(image1);
//     // x0.print_data();
//     // x1.print_data();

//     // image<type>::pointer im0_ptr(&image0);
//     // image<type>::pointer im1_ptr(&image1); 
//     grid<type>::pointer x0_ptr(&x0);
//     grid<type>::pointer x1_ptr(&x1);

//     // interpolate<type> image1_p(im1_ptr, x1_ptr); //*** CHANGE back interpolate to regular objects Mar 06 2020
//     // interpolate<type> image0_p(im0_ptr, x0_ptr);

//     interpolate<type> image0_p(image0, x0);
//     interpolate<type> image1_p(image1, x1);

//     // image<type> ssd_(image0.get_width(),image0.get_height());
//     // ssd_ = (image1_p*(transform * x0));
//     // image0 - image1_p*(transform * x0);
//     image<type> ssd_;
//     // ssd_ = image0 - image1_p*(transform * x0);
//     ssd_ = image1 - image0_p*(transform * x1);

//     return ssd_;
// };


image<type> ssd(image<type>::pointer image0, image<type>::pointer image1, affine<type> & transform)
{

    grid<type>::pointer x0 = grid<type>::new_pointer(*image0);
    grid<type>::pointer x1 = grid<type>::new_pointer(*image1);

    interpolate<type> image0_p(image0, x0);
    interpolate<type> image1_p(image1, x1);

    image<type> ssd_;
    // ssd_ = *image0 - image1_p*(transform * (*x0));
    ssd_ = *image1 - image0_p*(transform * (*x1));

    return ssd_;
};

int main()
{
    // ============================================
    //          Testing interpolate
    // ============================================

    // Create small imame
    // image<type>::pointer image0 = image<type>::new_pointer(4,3);
    image<type>::pointer image0 = image<type>::new_pointer(10,7);
    image0->random();
    // image0->ones();
    grid<type>::pointer x0 = grid<type>::new_pointer(*image0);

    // image<type> params(6,1);
    // *params.get_data() = {1.0, 0.0, 0.0, 1.0, 3.5, -1.0};
    image<type>::pointer params( new image<type>{1.1, 0.5, -0.5, 0.9, 2.1, -1.1} );
    affine<type> translation( 2, params);

    grid<type>::pointer t_x0 = grid<type>::new_pointer();
    *t_x0 = translation.transform(*x0);
    x0->print();
    t_x0->print();

    t_x0->print_data();
    translation.print();

    std::cout << "Create interpolation object\n";
    interpolate<type> image1_itp(image0, x0);
    std::cout << "Linear interpolation function\n";

    image<type>::pointer image0t = image<type>::new_pointer();
    *image0t = image1_itp.linear(*t_x0); 
    // TODO: Error of free() when using double. The error is in this method
    // *** TEST neighbors4??? change array to vector
    
    std::cout << "Linear interpolation finished\n";

    image0->print();
    image0->print_data();
    
    image0t->print();
    image0t->print_data();

    // Testing the inverse transfor applied to ssd function
    std::cout << "Test compressed notation\n";
    image<type>::pointer params_inv( new image<type>{1.0, 0.0, 0.0, 1.0, -3.5, 1.0});
    affine<type> inv_translation(2,params_inv);

    // image<type> result = ssd(image0, image0t, inv_translation);
    image<type> result = ssd(image0, image0t, translation);
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

    //     image<type> result = (image0 - image1_p*transform*x0)^2.0;
    //     return result.sum()
    // };

    return 0;

};