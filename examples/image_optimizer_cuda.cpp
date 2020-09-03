/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-29 13:11:53
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
    auto img = image_cuda<unsigned short>::new_pointer();
    img->read(filename);

    auto image0_cast = image_cuda<type>::new_pointer();
    auto image0 = image_cuda<type>::new_pointer();
    cast(*img, *image0_cast);

    // std::cout << "Min: " << image0_cast->min() << std::endl;
    // std::cout << "Max: " << image0_cast->max() << std::endl;

    *image0 = normalize(*image0_cast);
    // image0->print_data();

    // std::cout << "Min: " << image0->min() << std::endl;
    // std::cout << "Max: " << image0->max() << std::endl;

    auto x0 = grid<type,vector_cuda<type>>::new_pointer(*image0);

    image_cuda<type>::pointer params( new image_cuda<type>{1.0, 0.2, 0.0, 1.0, -20.0, 20.0} );
    auto taffine = affine<type,vector_cuda<type>>::new_pointer(2, params);
    // taffine->print_data();

    // auto imaget = image_cuda<type>::new_pointer();
    // auto image0_itpw = ilinear<type,vector_cuda<type>>::new_pointer(image0_cast);
    // imaget = image0_itpw->apply(taffine->apply(x0));
    // auto img_out = image_cuda<unsigned short>::new_pointer();
    // cast(*imaget, *img_out);
    // img_out->write("./transformed.png");

    auto image1 = image_cuda<type>::new_pointer();
    auto image0_itp = ilinear<type,vector_cuda<type>>::new_pointer(image0);
    image1 = image0_itp->apply(taffine->apply(x0));
    // image0->print();
    // image0->print_data();
    
    // image1->print(); 
    // image1->print_data();

    // Gradient
    // typename image_cuda<type>::vector grad(image0->get_dimension());
    // grad = gradient(image0);
    // image0->print_data("i");
    // grad[0]->print();
    // grad[0]->print_data("di/dx");
    // grad[1]->print_data("di/dy");

    auto trfm = affine<type,vector_cuda<type>>::new_pointer(2);
    // trfm->print_data();

    auto ssd1 = ssd<type,vector_cuda<type>>::new_pointer(image0, image1, trfm);
    // std::cout << ssd1->cost();

    // Affine test
    // *trfm + *taffine;
    // affine<type,vector_cuda<type>> tt = (*trfm - *taffine);
    // tt.print();
    // tt.print_data();

    auto opt = gradient_descent<type,vector_cuda<type>>::new_pointer();
    opt->optimize(ssd1);
    ssd1->get_transform()->print_data();
    return 0;
};