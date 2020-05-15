/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-04-18 13:47:14
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
#include "../src/optimizer.h"


// using type = int;
using type = float;
// using type = double;

int main()
{
    // ============================================
    //          Testing metric
    // ============================================

    // Create small imame
    std::string filename = "./tests/images/sinc_pad.png";
    auto img = image<unsigned short>::new_pointer();
    img->read(filename);

    auto image0_cast = image<type>::new_pointer();
    auto image0 = image<type>::new_pointer();
    type value = 0.0;
    *image0_cast = cast(*img, value);

    std::cout << "Min: " << image0_cast->min() << std::endl;
    std::cout << "Max: " << image0_cast->max() << std::endl;

    *image0_cast = normalize<type>(*image0_cast);
    image0->copy(*image0_cast);

    std::cout << "Min: " << image0->min() << std::endl;
    std::cout << "Max: " << image0->max() << std::endl;

    auto x0 = grid<type>::new_pointer(*image0);

    image<type>::pointer params( new image<type>{1.0, 0.2, 0.0, 1.0, -20.0, 20.0} );
    auto taffine = affine<type>::new_pointer(2, params);
    // taffine->print();

    auto image1 = image<type>::new_pointer();
    auto image0_itp = interpolate<type>::new_pointer(image0, x0);
    *image1 = image0_itp->linear(taffine->transform(*x0));

    unsigned short ww = 0;
    auto img_out = image<unsigned short>::new_pointer();
    *img_out = cast(*image1, ww);

    img_out->write("./transformed.png");

    image0->print();
    // image0->print_data();
    
    image1->print(); 
    // image1->print_data();

    typename image_base<type>::vector grad(image0->get_dimension());
    grad = gradient(*image0);
    // grad[0]->print_data();
    // grad[1]->print_data();

    image<type>::pointer iden( new image<type>{1.0, 0.0, 0.0, 1.0, 0.0, 0.0} );
    auto trfm = affine<type>::new_pointer(2, iden);
    // trfm->print_data();

    ssd<type> ssd1(image0, image1, trfm);
    // std::cout << ssd1.cost();
    optimizer<type> optim;
    optim.optimize(ssd1);
    
    return 0;

};