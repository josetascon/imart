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
    // Type definitions
    // using vector_image = std::vector<image_2d<pixel_type>>;

protected:
    
    // std::shared_ptr<vector_image> grid;
    // std::vector<image_2d> grid;

public:

    // std::shared_ptr<vector_image> get_grid();

    
    using image_base<pixel_type>::image_base;

    using image_base<pixel_type>::operator =;
    using image_base<pixel_type>::operator +;
    using image_base<pixel_type>::operator -;
    using image_base<pixel_type>::operator *;
    using image_base<pixel_type>::operator /;

    using image_base<pixel_type>::operator ();

    

    // TODO
    // assigment of index is not working e.g. image_2d(0,0) = 7
};


// ===========================================
//      Functions of Class image_2d
// ===========================================






#endif