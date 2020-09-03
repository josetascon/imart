/*
* @Author: Jose Tascon
* @Date:   2019-11-07 10:13:08
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-28 10:01:13
*/

// std libs
#include <iostream>
#include <vector>

// local libs
#include "../src/image.h"

using namespace imart;

int main()
{
    using type = double;
    
    // =====================================================
    // Two dimensional gradient
    // =====================================================
    // Create image
    int w = 4;
    int h = 5;
    auto img = image_cuda<type>::new_pointer(w,h);

    // Initialize
    std::vector<type> vec(w*h);
    type * p = vec.data();
    for(int i = 0; i < w*h; i++)
    {
        p[i] = i;
    };
    img->get_data()->read_ram(p,w*h);
    // img->random();
    
    // Comput gradient
    typename image_cuda<type>::vector grad(2);
    grad = gradient(img);
    img->print_data("i");
    grad[0]->print_data("didx");
    grad[1]->print_data("didy");

    // =====================================================
    // Three dimensional gradient
    // =====================================================
    // Create image
    int l = 3;
    auto img3 = image_cuda<type>::new_pointer(w,h,l);

    // Initialize
    std::vector<type> vec3(w*h*l);
    type * p3 = vec3.data();
    for(int i = 0; i < w*h*l; i++)
    {
        p3[i] = i;
    };
    img3->get_data()->read_ram(p3,w*h*l);
    // img3->random();
    
    // Comput gradient
    typename image_cuda<type>::vector grad3(3);
    grad3 = gradient(img3);
    img3->print_data("i");
    grad3[0]->print_data("didx");
    grad3[1]->print_data("didy");
    grad3[2]->print_data("didz");

    return 0;
}