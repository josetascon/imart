/*
* @Author: Jose Tascon
* @Date:   2020-06-06 00:00:00
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-08 10:14:47
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/vector_cpu.h"

using namespace imart;

using type = float;

int main()
{
    // ============================================
    //          Testing 
    // ============================================
    vector_cpu<type>::pointer vec1 = vector_cpu<type>::new_pointer();
    vector_cpu<type>::pointer vec2 = vector_cpu<type>::new_pointer(3);
    vector_cpu<type>::pointer vec3 = vector_cpu<type>::new_pointer(5, 1.1);
    vector_cpu<type>::pointer vec4 = vector_cpu<type>::new_pointer(*vec3);

    vec1->print();
    vec1->print_data();
    vec2->print_data();
    vec3->print_data();
    vec4->print_data();
    std::cout << std::endl;

    vector_cpu<type>::pointer vec11 = vector_cpu<type>::new_pointer(3,3.1);
    vector_cpu<type>::pointer vec12 = vector_cpu<type>::new_pointer(3,-0.2);

    vector_cpu<type>::pointer vec13 = vector_cpu<type>::new_pointer();
    vec13 = *vec11 + *vec12;

    vec11->print();
    vec11->print_data();
    vec12->print_data();
    vec13->print_data();

    // std::vector<type> vv1(3,1.2);
    // std::vector<type> vv2(3);
    // vv2 = vv1;

    // std::cout << "vector\n";
    // for(int k=0; k<3;k++){std::cout << vv2.at(k) << " ";};



    return 0;

};