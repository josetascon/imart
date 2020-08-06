/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-07-23 10:12:24
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/resolution.h"

using namespace imart;

int main()
{
    // ============================================
    //              Testing resolution
    // ============================================
    // using type = int;
    using type = float;
    // using type = double;
    
    // Create small imame
    auto image0 = image_cpu<type>::new_pointer(10,8);
    image0->set_spacing(std::vector<double>{2.0,1.5});
    image0->ones();
    
    auto mresolution1 = resolution<type,vector_cpu<type>>::new_pointer(image0);
    auto image1 = mresolution1->apply(4.0);

    image0->print();
    image0->print_data();
    image1->print();
    image1->print_data();

    // Create small imame
    auto image2 = image_gpu<type>::new_pointer(10,8);
    image2->set_spacing(std::vector<double>{2.0,1.5});
    image2->ones();
    
    auto mresolution2 = resolution<type,vector_ocl<type>>::new_pointer(image2);
    auto image3 = mresolution2->apply(4.0);

    image2->print();
    image2->print_data();
    image3->print();
    image3->print_data();

    return 0;
};