/*
* @Author: Jose Tascon
* @Date:   2020-06-19 20:05:07
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-23 05:30:01
*/

// std libs
#include <iostream>

// local libs
#include "../src/opencl_object.h"
#include "../src/vector_opencl.h"

using namespace imart;

int main()
{
    using type = double;

    vector_opencl<type> veca;
    veca.print();
    veca.print_data();

    vector_opencl<type> vecb(10);
    vecb.print();
    vecb.print_data();

    int N = 20;
    type scalar = 3;
    vector_opencl<type> vecc(N,scalar);
    vecc.print();
    vecc.print_data();

    auto vec1 = vector_opencl<type>::new_pointer(N,4);
    auto vec2 = vec1->clone();
    auto vec3 = vec1->copy();
    auto vec4 = vec1->mimic();

    vec1->print("vec1");
    vec1->print_data();
    vec2->print("vec2 = vec1->clone()");
    vec2->print_data();
    vec3->print("vec3 = vec1->copy()");
    vec3->print_data();
    vec4->print("vec4 = vec1->mimic()");
    vec4->print_data();

    auto vec11 = vector_opencl<type>::new_pointer(10,3.1);
    auto vec12 = vector_opencl<type>::new_pointer(10,-0.2);

    auto vec13 = vector_opencl<type>::new_pointer();
    vec13 = *vec11 + *vec12;

    vec11->print();
    vec11->print_data();
    vec12->print_data();
    vec13->print_data();

    std::cout << "vec12[0]: ";
    std::cout << vec12->operator[](0) << std::endl;
    std::cout << "vec12[9]: ";
    std::cout << vec12->operator[](9) << std::endl;

    std::vector<type> std_vec(22);
    for(int i = 0; i < std_vec.size(); i++)
    {
        std_vec[i] = i*2.2;
    }

    std::cout << "std_vec:\n[ ";
    for(int i = 0; i < std_vec.size(); i++) std::cout << std_vec[i] << " ";
    std::cout << "]" << std::endl;

    auto vec_opencl = vector_opencl<type>::new_pointer(22);
    vec_opencl->print_data("create vec_opencl");
    vec_opencl->read_ram(std_vec.data(), std_vec.size());
    vec_opencl->print_data("copy ram std_vec");

    int ss = 501;
    std::vector<type> std_vec1(ss);
    type add = 0;
    for(int i = 0; i < std_vec1.size(); i++)
    {
        std_vec1[i] = 101.2 - i;
        add += 101.2 - i;
    }

    auto vec_opencl1 = vector_opencl<type>::new_pointer(ss);
    vec_opencl1->read_ram(std_vec1.data(), std_vec1.size());
    // vec_opencl1->print_data();
    std::cout << "add: " << add << std::endl;
    std::cout << "add opencl: " << vec_opencl1->sum() << std::endl;
    std::cout << "min opencl: " << vec_opencl1->min() << std::endl;
    std::cout << "max opencl: " << vec_opencl1->max() << std::endl;


    return 0;
}

    