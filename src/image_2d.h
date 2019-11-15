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
#include <assert.h>       // assert

// local lib
#include "image_base.h"



// Class image_base
template <typename pixel_type>
class image_2d : public image_base<pixel_type>
{
public:
    using image_base<pixel_type>::image_base;

    using image_base<pixel_type>::operator=;
    using image_base<pixel_type>::operator+;
    using image_base<pixel_type>::operator-;
    using image_base<pixel_type>::operator*;
    using image_base<pixel_type>::operator/;
};



#endif