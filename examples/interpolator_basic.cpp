/*
* @Author: Jose Tascon
* @Date:   2020-06-23 12:30:38
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-23 12:43:54
*/

// std libs
#include <iostream>     // std::cout

// local libs
#include "../src/interpolator.h"

using namespace imart;

int main()
{
    using type = float;

    auto img1 = image<type>::new_pointer(5,3);
    auto x1 = grid<type>::new_pointer(*img1);

    auto itp = interpolator<type>::new_pointer(img1, x1);
    itp->print();

    return 0;
};