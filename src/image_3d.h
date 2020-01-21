/*
* @Author: jose
* @Date:   2019-11-26 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-26 00:00:00
*/

#ifndef __IMAGE_3D_H__
#define __IMAGE_3D_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector

// local lib
#include "image_base.h"


// Class image_base
template <typename pixel_type>
class image_3d : public image_base<pixel_type>
{
public:
    // Type definitions
    using ptr_pixels8 = std::unique_ptr<std::array<pixel_type,8>>;
    
    image_3d();
    image_3d(int w, int h, int l);
    image_3d(std::shared_ptr<std::vector<pixel_type>> buffer, int w, int h, int l) : image_base<pixel_type> (buffer,w,h,l) { ; };
    image_3d(const image_3d<pixel_type> & input) : image_base<pixel_type> (input) { ; };

    using image_base<pixel_type>::operator =;

    // std::string info_data(std::string msg); //ONLY 3d***


    // Access
    ptr_pixels8 neighbors8(int e);

    // TODO
    // index assignment is not working image_3d(0,0,0) = 7
};


// ===========================================
//      Functions of Class image_3d
// ===========================================
template <typename pixel_type>
image_3d<pixel_type>::image_3d()
{
    image_base<pixel_type>::init(0, 0, 0);
    this->class_name = "image_3d";
    this->data.reset();
};

template <typename pixel_type>
image_3d<pixel_type>::image_3d(int w, int h, int l): image_base<pixel_type> (w,h,l)
{
    this->class_name = "image_3d";
};

// template <typename pixel_type>
// std::string image_3d<pixel_type>::info_data(std::string msg)
// {
//     std::stringstream ss;
//     if (msg != "") { ss << msg << std::endl; };
//     int w = this->width;
//     int h = this->height;
//     int l = this->length;
//     pixel_type * p = this->ptr();

//     // std::cout << "Image data:" << std::endl;
//     // ss << "[";
//     for(int k = 0; k < l; k++)
//     {
//         ss << "[ ";
//         for(int i = 0; i < h; i++)
//         {
//             for(int j=0; j < w; j++)
//             {
//                 ss << p[j + i*w + k*w*h] << " "; // valgrind error solved
//             };
//             if(i < h-1){ss << std::endl << "  ";};
//         };
//         ss << "]" << std::endl;
//     };
//     // ss << "]";
//     ss << std::endl;
//     return ss.str();
// };


template <typename pixel_type>
typename image_3d<pixel_type>::ptr_pixels8 image_3d<pixel_type>::neighbors8(int e)
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