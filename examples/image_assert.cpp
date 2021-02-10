/*
* @Author: jose
* @Date:   2019-11-05 14:55:42
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-01 22:59:13
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

    // int N = 4*3*2+1;

    // for(int i = 0; i < N; i++)
    // {
    //     std::cout << img1->operator[](i) << " ";
    // };


    // Test equal function
    auto img4 = image<type>::new_pointer(4,3,2);
    auto img5 = image<type>::new_pointer(4,3,2);

    img5->print(); // check pointer befor equal operation

    img4->random();
    img5->equal(*img4);

    img4->print();
    img4->print_data();
    img5->print();
    img5->print_data();



    return 0;
};