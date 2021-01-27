/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-10-08 16:11:36
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"

using namespace imart;

int main()
{
    using type = float;
    auto img1 = image<type>::new_pointer(4,3,2);
    auto img2 = image<type>::new_pointer(5,3,2);
    auto img3 = image<type>::new_pointer(3);

    img1->ones();
    img2->ones();

    // *img3 = *img1 + *img2;

    int N = 4*3*2+1;

    for(int i = 0; i < N; i++)
    {
        std::cout << img1->operator[](i) << " ";
    };


    return 0;
};