/*
* @Author: Jose Tascon
* @Date:   2020-06-19 20:05:07
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-21 14:48:46
*/

// std libs
#include <iostream>

// local libs
#include "../src/image.h"
#include "../src/image_utils.h"

using namespace imart;

int main()
{
    using type = double;

    auto imgc = image_cpu<type>::new_pointer(4,5);
    auto imgg = image_gpu<type>::new_pointer(4,5);
    std::vector<type> vec{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    
    imgc->get_data()->read_ram(vec.data(),vec.size());
    imgg->get_data()->read_ram(vec.data(),vec.size());

    auto gxc = gradientx(imgc);
    auto gyc = gradienty(imgc);
    imgc->print_data("image cpu");
    gxc->print_data("gx");
    gyc->print_data("gy");

    auto gxg = gradientx(imgg);
    auto gyg = gradienty(imgg);
    imgg->print_data("image gpu");
    gxg->print_data("gx");
    gyg->print_data("gy");


    return 0;
}

    