/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-06 19:21:51
*/

// std libs
#include <iostream>
#include <memory>
#include <vector>

// local libs
#include "../src/image.h"
#include "../src/affine.h"

using namespace imart;

int main()
{
    using type = float;

    // ============================================
    //              Testing 2d affine
    // ============================================
    // Create transform
    image_cpu<type>::pointer params1(new image_cpu<type>{1.1, 0.5, -0.5, 0.9, 2.1, -1.1});
    image_cpu<type>::pointer params2(new image_cpu<type>{-0.1 + 1.0, -0.5, 0.5, 1.1, -2.1 + 10.0, 1.1 - 10.0});
    auto affine1 = affine_cpu<type>::new_pointer(2,params1);
    auto affine2 = affine_cpu<type>::new_pointer(2,params2);
    auto affine3 = affine_cpu<type>::new_pointer(2);
    auto affine4 = affine_cpu<type>::new_pointer(2);

    *affine3 = *affine1 + *affine2;
    
    affine1->print_data("affine1:");
    affine2->print_data("affine2:");
    affine3->print_data("affine3 = affine1 + affine2");

    *affine4 = (*affine3)*0.5;

    affine4->print_data("affine4 = affine3*0.5");

    image_cpu<type>::pointer params3(new image_cpu<type>{3.0, 3.0, 3.0, 3.0, 3.0, 3.0});
    auto affine5 = affine_cpu<type>::new_pointer(2, params3);
    auto affine6 = affine_cpu<type>::new_pointer(2);

    *affine6 = (*affine4)*(*affine5);

    affine5->print_data("affine5:");
    affine6->print_data("affine6 = affine4*affine5");

    

    return 0;
};