/*
* @Author: jose
* @Date:   2019-11-05 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __IMAGE_H__
#define __IMAGE_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector

// local lib
#include "image_base.h"
// #include "utils/timer.h"


// Class image_base
template <typename pixel_type>
class image : public image_base<pixel_type>
{
public:
    //Type definitions
    using self    = image;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<image::pointer>;

protected:
    // Type definitions
    using ptr_pixels4 = std::unique_ptr<std::array<pixel_type,4>>;
    using ptr_pixels8 = std::unique_ptr<std::array<pixel_type,8>>;

public:

    image() : image_base<pixel_type> () { this->class_name = "image"; };
    image(int d) : image_base<pixel_type> (d) { this->class_name = "image"; };
    image(int w, int h) : image_base<pixel_type> (w,h) { this->class_name = "image"; };
    image(int w, int h, int l) : image_base<pixel_type> (w,h,l) { this->class_name = "image"; };
    image(std::shared_ptr<std::vector<pixel_type>> buffer, int w, int h) : image_base<pixel_type> (buffer,w,h) { this->class_name = "image"; };
    image(std::shared_ptr<std::vector<pixel_type>> buffer, int w, int h, int l) : image_base<pixel_type> (buffer,w,h,l) { this->class_name = "image"; };
    image(std::initializer_list<pixel_type> list) : image_base<pixel_type> (list) { this->class_name = "image"; };
    image(const image<pixel_type> & input) : image_base<pixel_type> (input) { this->class_name = "image"; };
    image(const image_base<pixel_type> & input) : image_base<pixel_type> (input) { this->class_name = "image"; };

    using image_base<pixel_type>::operator =;
    // using image_base<pixel_type>::operator +;
    // using image_base<pixel_type>::operator -;
    
    // Access
    ptr_pixels4 neighbors4(int e);
    ptr_pixels8 neighbors8(int e);

    // TODO
    // assigment of index is not working e.g. image_2d(0,0) = 7
};


// ===========================================
//      Functions of Class image
// ===========================================
template <typename pixel_type>
typename image<pixel_type>::ptr_pixels4 image<pixel_type>::neighbors4(int e)
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

template <typename pixel_type>
typename image<pixel_type>::ptr_pixels8 image<pixel_type>::neighbors8(int e)
{
    ptr_pixels8 arr = std::make_unique<std::array<pixel_type,8>>();
    pixel_type * p = this->ptr();
    int w = this->width;
    int h = this->height;

    (*arr)[0] = p[e];
    (*arr)[1] = p[e+1];
    (*arr)[2] = p[e+w];
    (*arr)[3] = p[e+w+1];
    
    (*arr)[4] = p[e + w*h];
    (*arr)[5] = p[e+1 + w*h];
    (*arr)[6] = p[e+w + w*h];
    (*arr)[7] = p[e+w+1 + w*h];
    // std::cout << "px:" << p[e] << " ";
    return arr;
};



#endif