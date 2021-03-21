/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-07-28 13:01:39
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/ssd.h"
#include "../src/cc.h"

using namespace imart;

int main()
{
    // ============================================
    //              Testing ssd metric
    // ============================================
    // using type = int;
    using type = float;
    // using type = double;
    
    // Create small imame
    auto image0 = image<type,vector_opencl<type>>::new_pointer(10,7);
    auto image1 = image0->mimic();
    image0->assign(1.0);
    image1->assign(2.0);

    auto ssd1 = ssd<type,vector_opencl<type>>::new_pointer(image0, image1);
    // auto ssd1 = ssd<type>::new_pointer();
    ssd1->print();
    image0->print();
    image0->print_data("\nfixed image");
    image1->print_data("moving image");
    std::cout << "ssd cost: " << ssd1->cost() << std::endl;
    // std::cout << "ssd cost: " << ssd1->cost(image0, image1) << std::endl;

    // ============================================
    //              Testing cc metric
    // ============================================
    image0->random();
    image1->random();
    auto cc1 = cc<type,vector_opencl<type>>::new_pointer(image0, image1);
    std::cout << "cc cost: " << cc1->cost() << std::endl;
    
    return 0;
};