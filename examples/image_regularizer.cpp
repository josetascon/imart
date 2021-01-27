/*
* @Author: jose
* @Date:   2021-01-20 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-21 19:06:01
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/regularizer.h"

using namespace imart;

int main()
{
    using type = float;
    using container = vector_cpu<type>;
    
    auto reg = regularizer<type,container>::new_pointer(1.0,1.0);

    // Test K
    auto v = std::make_shared< std::vector<std::shared_ptr<image<type,container>>> >(2);
    
    auto i0 = image<type,container>::new_pointer(5,5);
    auto i1 = image<type,container>::new_pointer(5,5);
    
    auto dat = vector_cpu<type>::new_pointer(5*5, 0.0);

    dat->data()[2+2*5] = 1.0;
    i0->set_data(dat->clone());

    dat->data()[2+2*5] = 0.0;
    i1->set_data(dat->clone());
    
    v->at(0) = i0;
    v->at(1) = i1;

    reg->update_a(i0);
    auto gv = reg->k(v);
    std::cout << "Regularizer K" << std::endl;
    std::cout << "Input" << std::endl;
    i0->print_data();
    i1->print_data();
    std::cout << "K" << std::endl;
    gv[0]->print_data();
    gv[1]->print_data();
    

    // Test A
    int w = 5;
    int h = 4;

    auto img = image<type,container>::new_pointer(w,h);

    auto g = std::make_shared< std::vector<std::shared_ptr<image<type,container>>> >(2);
    g->at(0) = img;

    auto out = reg->update_a(img);

    std::cout << "Regularizer A" << std::endl;
    // out[0]->print();
    out[0]->print_data();
    out[1]->print_data();

    return 0;
};