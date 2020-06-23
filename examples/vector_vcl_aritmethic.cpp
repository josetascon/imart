/*
* @Author: Jose Tascon
* @Date:   2020-06-06 00:00:00
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-22 11:48:30
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/vector_vcl.h"

using namespace imart;

using type = float;

int main()
{
    // ============================================
    //          Testing 
    // ============================================
    vector_vcl<type>::pointer vec1 = vector_vcl<type>::new_pointer();
    vector_vcl<type>::pointer vec2 = vector_vcl<type>::new_pointer(3);
    vector_vcl<type>::pointer vec3 = vector_vcl<type>::new_pointer(10, 2.1);
    vector_vcl<type>::pointer vec4 = vector_vcl<type>::new_pointer(*vec3);
    vector_vcl<type>::pointer vec5 = vec4->clone();
    vector_vcl<type>::pointer vec6 = vector_vcl<type>::new_pointer();

    vec1->print();
    vec1->print_data("vec1:");
    vec2->print_data("vec2:");
    vec3->print_data("vec3:");
    vec4->print("vec4");
    vec4->print_data("vec4:");
    vec5->print("vec5");
    vec4->print_data("vec5:");
    vec6 = *vec5 + *vec5;
    vec6->print_data("vec6");
    std::cout << std::endl;

    vector_vcl<type>::pointer vec11 = vector_vcl<type>::new_pointer(10,3.1);
    vector_vcl<type>::pointer vec12 = vector_vcl<type>::new_pointer(10,-0.2);

    vector_vcl<type>::pointer vec13 = vector_vcl<type>::new_pointer();
    vec13 = *vec11 + *vec12;

    vec11->print();
    vec11->print_data();
    vec12->print_data();
    vec13->print_data();

    int s = 20;
    vector_vcl<type>::pointer vec21 = vector_vcl<type>::new_pointer(s,1.0);
    vector_vcl<type>::pointer vec22 = vector_vcl<type>::new_pointer(s,1.0);
    vector_vcl<type>::pointer vec25 = vector_vcl<type>::new_pointer(s,1.0);
    vector_vcl<type>::pointer vec23 = vector_vcl<type>::new_pointer();
    vector_vcl<type>::pointer vec24 = vector_vcl<type>::new_pointer();
    vector_vcl<type>::pointer vec26 = vector_vcl<type>::new_pointer();

    std::vector<type> a(s); 
    std::vector<type> b(s);
    std::vector<type> c(s);
    for (int i=0; i<s; ++i)
    {
        // vec21->operator[](i) = static_cast<type>(-i-131.2);
        // vec22->operator[](i) = static_cast<type>(i+10.5);
        a[i] = type(-i-131.2);
        b[i] = type(i+10.5);
        c[i] = type(i);
    };

    viennacl::copy(a.begin(), a.end(), vec21->begin());
    viennacl::copy(b.begin(), b.end(), vec22->begin());
    viennacl::copy(c.begin(), c.end(), vec25->begin());

    vec23 = vec21->normalize();
    vec23->print();
    vec24 = vec22->normalize();
    vec24->print();
    vec26 = vec24->normalize();
    vec26->print();

    // vec23 = *vec21-2.0;
    // vec23->print();
    // vec24 = *vec21-3.0;
    // vec24->print();
    // vec26 = *vec21-4.0;
    // vec26->print();

    vec23->print_data();
    vec24->print_data();
    vec26->print_data();


    return 0;

};