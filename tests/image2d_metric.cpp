/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-03-19 15:48:46
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"
#include "../src/interpolate.h"
#include "../src/ssd.h"


// using type = int;
using type = float;
// using type = double;

int main()
{
    // ============================================
    //          Testing metric
    // ============================================

    // Create small imame
    auto image0 = image<type>::new_pointer(10,7);
    image0->random();
    // image0.ones();
    auto x0 = grid<type>::new_pointer(*image0);

    image<type>::pointer params( new image<type>{1.1, 0.5, -0.5, 0.9, 2.1, -1.1} );
    // affine<type> taffine( 2, params);
    // taffine.print();
    auto taffine = affine<type>::new_pointer(2, params);
    taffine->print();

    auto image0t = image<type>::new_pointer();
    auto image0_itp = interpolate<type>::new_pointer(*image0, *x0);
    *image0t = image0_itp->linear(taffine->transform(*x0));

    image0->print();
    image0->print_data();
    // std::cout << image0->ptr() << std::endl;
    
    image0t->print();
    image0t->print_data();
    // std::cout << image0t->ptr() << std::endl;

    
    ssd<type> ssd1(image0, image0t, taffine);
    std::cout << "SSD: " << ssd1.cost() << std::endl;

    ssd1.print();
    ssd1.print_data();

    // Old style, everything by reference
    // // Create small imame
    // image<type> image0(10,7);
    // image0.random();
    // // image0.ones();
    // grid<type> x0(image0);

    // image<type>::pointer params( new image<type>{1.1, 0.5, -0.5, 0.9, 2.1, -1.1} );
    // // affine<type> taffine( 2, params);
    // // taffine.print();
    // affine<type>::pointer taffine(new affine<type>(2, params));
    // taffine->print();

    // interpolate<type> image0_itp(image0, x0);
    // image<type> image0t = image0_itp*((*taffine)*x0);

    // image0.print();
    // image0.print_data();
    
    // image0t.print();
    // image0t.print_data();

    
    // ssd<type> ssd1(image0, image0t, taffine);
    // std::cout << ssd1.cost() << std::endl;

    // ssd1.print();
    // ssd1.print_data();

    return 0;

};