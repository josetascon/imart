/*
* @Author: jose
* @Date:   2019-11-19 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-19 13:55:00
*/

#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <cassert>      // assert

// images 
#include "image_base.h"
#include "grid.h"


// Class image_base_2d
template <typename pixel_type>
class interpolate
{
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int dim;
    std::string type;
    image_base<pixel_type> image;
    grid<pixel_type> xyz;
    
public:
    // ===========================================
    // Create Functions
    // ===========================================
    // interpolate();
    interpolate(image_base<pixel_type> & input, grid<pixel_type> xref);

    // ===========================================
    // Functions
    // ===========================================
    image_base<pixel_type> linear(grid<pixel_type> xin);
};


template <typename pixel_type>
interpolate<pixel_type>::interpolate(image_base<pixel_type> & input, grid<pixel_type> xref)
{
    image = input;
    
};

