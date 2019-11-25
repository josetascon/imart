/*
* @Author: jose
* @Date:   2019-11-05 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __IMAGE_2D_H__
#define __IMAGE_2D_H__

// std libs
#include <iostream>     // std::cout
#include <sstream>      // stringstream
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <random>       // random
#include <typeinfo>     // operator typeid
#include <cassert>      // assert

// local lib
#include "image_base.h"
// #include "utils/timer.h"


// Class image_base
template <typename pixel_type>
class image_2d : public image_base<pixel_type>
{
public:
    // Type definitions
    // using vector_image = std::vector<image_2d<pixel_type>>;
    using ptr_pixels4 = std::unique_ptr<std::array<pixel_type,4>>;
// protected:
    
    // std::shared_ptr<vector_image> grid;
    // std::vector<image_2d> grid;

public:
    
    image_2d() : image_base<pixel_type> () { ; };
    image_2d(int w, int h) : image_base<pixel_type> (w,h) { ; };
    image_2d(std::shared_ptr<std::vector<pixel_type>> buffer, int w, int h) : image_base<pixel_type> (buffer,w,h) { ; };
    image_2d(const image_2d<pixel_type> & input) : image_base<pixel_type> (input) { ; };

    // std::shared_ptr<vector_image> get_grid();

    // write this
    // using image_base<pixel_type>::image_base;

    // image_2d<pixel_type> & operator = (image_2d<pixel_type> input);
    // image_2d<pixel_type> operator + (image_2d<pixel_type> & input);
    using image_base<pixel_type>::operator =;

    // using image_base<pixel_type>::operator +;
    // using image_base<pixel_type>::operator -;
    // using image_base<pixel_type>::operator *;
    // using image_base<pixel_type>::operator /;

    // using image_base<pixel_type>::operator ();

    // using image_base<pixel_type>::print_data;

    // Access
    ptr_pixels4 neighbors4(int e);

    // TODO
    // assigment of index is not working e.g. image_2d(0,0) = 7
};


// ===========================================
//      Functions of Class image_2d
// ===========================================
template <typename pixel_type>
typename image_2d<pixel_type>::ptr_pixels4 image_2d<pixel_type>::neighbors4(int e)
{
    ptr_pixels4 arr = std::make_unique<std::array<pixel_type,4>>();
    pixel_type * p = this->ptr();
    (*arr)[0] = p[e];
    (*arr)[1] = p[e+1];
    (*arr)[2] = p[e+this->width];
    (*arr)[3] = p[e+this->width+1];
    // std::cout << "px:" << p[e] << " ";
    return arr;
};


// template <typename pixel_type>
// image_2d<pixel_type> & image_2d<pixel_type>::operator = (image_2d<pixel_type> input)
// {
//     // delete &data;
//     this->data.reset();
//     this->data = input.get_data();
//     this->update(input);
//     return *this;
// };

// template <typename pixel_type>
// image_2d<pixel_type> & image_2d<pixel_type>::operator = (image_2d<pixel_type> input)
// {
//     return image_base<pixel_type>::operator=(static_cast<image_base<pixel_type>>(input));
// };

// template <typename pixel_type>
// image_2d<pixel_type> image_2d<pixel_type>::operator + (image_2d<pixel_type> & input)
// {
//     return image_base<pixel_type>::operator+((input));
// };



#endif