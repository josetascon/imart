/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-02 08:11:14
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/affine.h"
#include "../src/ilinear.h"
#include "../src/ssd.h"
#include "../src/gradient_descent.h"

using namespace imart;

int main()
{
    using type = float;
    // using type = double;

    // ============================================
    //          Testing metric
    // ============================================

    // Create small imame
    std::string filename = "./examples/images/sinc_pad.png";
    auto img = image<unsigned short>::new_pointer();
    img->read(filename);

    auto image0_cast = image<type>::new_pointer();
    auto image0 = image<type>::new_pointer();
    cast(*img, *image0_cast);

    // std::cout << "Min: " << image0_cast->min() << std::endl;
    // std::cout << "Max: " << image0_cast->max() << std::endl;

    *image0 = normalize<type>(*image0_cast);

    // std::cout << "Min: " << image0->min() << std::endl;
    // std::cout << "Max: " << image0->max() << std::endl;

    auto x0 = grid<type>::new_pointer(*image0);

    image<type>::pointer params( new image<type>{1.0, 0.2, 0.0, 1.0, -20.0, 20.0} );
    auto taffine = affine<type>::new_pointer(2, params);
    // taffine->print_data();

    // auto imaget = image<type>::new_pointer();
    // auto image0_itpw = ilinear<type>::new_pointer(image0_cast);
    // imaget = image0_itpw->apply(taffine->apply(x0));
    // auto img_out = image<unsigned short>::new_pointer();
    // cast(*imaget, *img_out);
    // img_out->write("./transformed.png");

    auto image1 = image<type>::new_pointer();
    auto image0_itp = ilinear<type>::new_pointer(image0);
    image1 = image0_itp->apply(taffine->apply(x0));
    // image0->print();
    // image0->print_data();
    
    // image1->print(); 
    // image1->print_data();

    // Gradient
    // typename image<type>::vector grad(image0->get_dimension());
    // grad = gradient(image0);
    // // grad[0]->print_data();
    // // grad[1]->print_data();

    auto trfm = affine<type>::new_pointer(2);
    // trfm->print_data();

    auto ssd1 = ssd<type>::new_pointer(image0, image1, trfm);
    // std::cout << ssd1->cost();

    // Affine test
    // *trfm + *taffine;
    // affine<type> tt = (*trfm - *taffine);
    // tt.print();
    // tt.print_data();

    auto opt = gradient_descent<type>::new_pointer();
    opt->optimize(ssd1);
    return 0;
};