/*
* @Author: jose
* @Date:   2019-11-19 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-19 13:55:00
*/

#ifndef __INEAREST_H__
#define __INEAREST_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector

// images 
#include "interpolator.h"

namespace imart
{

// Class linear interpolator
template <typename type, typename container>
class inearest: public inherit<inearest<type,container>, interpolator<type,container>>
{
public:
    //Type definitions
    using self    = inearest;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    using interpolator<type,container>::_default_;
    using interpolator<type,container>::region;
    using interpolator<type,container>::image_reference;
    using interpolator<type,container>::image_output;
    using interpolator<type,container>::x_reference;
    using interpolator<type,container>::x_output;

    using inherit<inearest<type,container>, 
                  interpolator<type,container>>::inherit;

protected:
    
    // ===========================================
    // Functions
    // ===========================================
    void boundaries2();
    void boundaries3();
    int search_(type x, type y);
    int search_(type x, type y, type z);

    typename image<type,container>::pointer nearest2(const typename grid<type,container>::pointer xout);
    typename image<type,container>::pointer nearest3(const typename grid<type,container>::pointer xout);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    inearest(typename image<type,container>::pointer imgref, typename grid<type,container>::pointer xref);

    // ===========================================
    // Functions
    // ===========================================
    typename image<type,container>::pointer apply(const typename grid<type,container>::pointer xout);
};

template<typename type>
using inearest_cpu = inearest<type,vector_cpu<type>>;

template<typename type>
using inearest_gpu = inearest<type,vector_ocl<type>>;


// ===========================================
//      Functions of Class grid
// ===========================================

// ===========================================
// Constructor
// ===========================================
template <typename type, typename container>
inearest<type,container>::inearest(typename image<type,container>::pointer imgref, typename grid<type,container>::pointer xref)
{
    assert(imgref->get_dimension() == xref->get_dimension());
    
    this->class_name = "nearest interpolator";
    int d = imgref->get_dimension();
    init(d);
    image_reference = imgref;
    x_reference = xref;
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
typename image<type,container>::pointer inearest<type,container>::apply(const typename grid<type,container>::pointer xout)
{
    auto a = image<type,container>::new_pointer(xout->get_dimension());
    if(xout->get_dimension() == 2){ a = nearest2(xout); };
    if(xout->get_dimension() == 3){ a = nearest3(xout); };
    return a;
};


template <typename type, typename container>
typename image<type,container>::pointer inearest<type,container>::nearest2(const typename grid<type,container>::pointer xout)
{
    assert(xout->get_dimension() == image_reference->get_dimension());

    int w = image_reference->get_size()[0];
    int h = image_reference->get_size()[1];

    auto image_out = image<type,container>::new_pointer(w, h);
    image_out->set_sod_parameters(xout->get_spacing(), xout->get_origin(), xout->get_direction());
    image_out->assign(_default_);

    typename image<type,container>::pointer xo = (xout->ptr()[0]);   // raw pointer to
    typename image<type,container>::pointer yo = (xout->ptr()[1]);   // input grid coordinates

    // Future step
    // convert xout to x_reference

    container::nearest2(xo->get_data(),yo->get_data(),image_reference->get_data(),
                image_out->get_data(), image_reference->get_size(), xout->get_size() );
    // image_out->set_data(imgo);

    return image_out;
};

template <typename type, typename container>
typename image<type,container>::pointer inearest<type,container>::nearest3(const typename grid<type,container>::pointer xout)
{
    assert(xout->get_dimension() == image_reference->get_dimension());

    int w = image_reference->get_size()[0];
    int h = image_reference->get_size()[1];
    int l = image_reference->get_size()[2];

    auto image_out = image<type,container>::new_pointer(w, h, l);
    image_out->set_sod_parameters(xout->get_spacing(), xout->get_origin(), xout->get_direction());
    image_out->assign(_default_);

    typename image<type,container>::pointer xo = (xout->ptr()[0]);   // raw pointer to
    typename image<type,container>::pointer yo = (xout->ptr()[1]);   // input grid coordinates
    typename image<type,container>::pointer zo = (xout->ptr()[2]);

    // Future step
    // convert xout to x_reference

    container::nearest3(xo->get_data(),yo->get_data(), zo->get_data(), image_reference->get_data(),
                image_out->get_data(), image_reference->get_size(), xout->get_size() );
    // image_out->set_data(imgo);

    return image_out;
};

}; //end namespace

#endif