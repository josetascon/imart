/*
* @Author: jose
* @Date:   2019-11-19 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-19 13:55:00
*/

#ifndef __ILINEAR_H__
#define __ILINEAR_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector

// images 
#include "interpolator.h"

namespace imart
{

// Class linear interpolator
template <typename type, typename container=vector_cpu<type>>
class ilinear: public inherit<ilinear<type,container>, interpolator<type,container>>
{
public:
    //Type definitions
    using self    = ilinear;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using interpolator<type,container>::_default_;
    using interpolator<type,container>::image_reference;
    using interpolator<type,container>::image_output;
    using interpolator<type,container>::x_output;

    using inherit<ilinear<type,container>, 
                  interpolator<type,container>>::inherit;

protected:
    
    // ===========================================
    // Functions
    // ===========================================
    typename image<type,container>::pointer linear2(const typename grid<type,container>::pointer xout);
    typename image<type,container>::pointer linear3(const typename grid<type,container>::pointer xout);
    typename image<type,container>::pointer linear2_test(const typename grid<type,container>::pointer xout);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    ilinear() : inherit<ilinear<type,container>, interpolator<type,container>>()
                { this->class_name = "linear interpolator"; };
    ilinear(int d) : inherit<ilinear<type,container>, interpolator<type,container>>(d)
                { this->class_name = "linear interpolator"; };
    ilinear(typename image<type,container>::pointer imgref);

    // ===========================================
    // Functions
    // ===========================================
    typename image<type,container>::pointer apply(const typename grid<type,container>::pointer xout);
    typename image<type,container>::pointer apply_test(const typename grid<type,container>::pointer xout);
};

template<typename type>
using ilinear_cpu = ilinear<type,vector_cpu<type>>;

// template<typename type>
// using ilinear_gpu = ilinear<type,vector_ocl<type>>;

#ifdef IMART_WITH_OPENCL
template<typename type>
using ilinear_ocl = ilinear<type,vector_ocl<type>>;
#endif

#ifdef IMART_WITH_CUDA
template<typename type>
using ilinear_cuda = ilinear<type,vector_cuda<type>>;
#endif

// ===========================================
//      Functions of Class grid
// ===========================================

// ===========================================
// Constructor
// ===========================================
// template <typename type, typename container>
// ilinear<type,container>::ilinear(typename image<type,container>::pointer imgref, typename grid<type,container>::pointer xref)
template <typename type, typename container>
ilinear<type,container>::ilinear(typename image<type,container>::pointer imgref)
{
    // assert(imgref->get_dimension() == xref->get_dimension());
    this->class_name = "linear interpolator";
    int d = imgref->get_dimension();
    this->init(d);
    image_reference = imgref;
    // x_reference = xref;
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
typename image<type,container>::pointer ilinear<type,container>::apply(const typename grid<type,container>::pointer xout)
{
    auto a = image<type,container>::new_pointer(xout->get_dimension());
    if(xout->get_dimension() == 2){ a = linear2(xout); }
    else if(xout->get_dimension() == 3){ a = linear3(xout); };
    // update process variables
    // x_output = xout; // update x_output
    // image_out = a;  // update image_out
    return a;
};

template <typename type, typename container>
typename image<type,container>::pointer ilinear<type,container>::apply_test(const typename grid<type,container>::pointer xout)
{
    auto a = image<type,container>::new_pointer(xout->get_dimension());
    if(xout->get_dimension() == 2){ a = linear2_test(xout); }
    return a;
};

template <typename type, typename container>
typename image<type,container>::pointer ilinear<type,container>::linear2_test(const typename grid<type,container>::pointer xout)
{
    assert(xout->get_dimension() == image_reference->get_dimension());

    int w = xout->get_size()[0];
    int h = xout->get_size()[1];

    // setup output image
    auto image_out = image<type,container>::new_pointer(w, h);
    image_out->set_sod_parameters(xout->get_spacing(), xout->get_origin(), xout->get_direction());
    // image_out->assign(_default_);

    typename image<type,container>::pointer xo = (xout->ptr()[0]);   // raw pointer to
    typename image<type,container>::pointer yo = (xout->ptr()[1]);   // input grid coordinates

    // convert xout to x_reference (world to image coordinates)
    // std::cout << "convert xout to xref" << std::endl;
    std::vector<double> sod = image_reference->get_sod_inverse();
    
    // apply the interpolation algorithm (or kernel)
    // std::cout << "linear call to container" << std::endl;
    container::linear2_test(xo->get_data(),yo->get_data(),image_reference->get_data(),
                image_out->get_data(), image_reference->get_size(), sod, _default_ );
    // image_out->set_data(imgo);

    return image_out;
};


template <typename type, typename container>
typename image<type,container>::pointer ilinear<type,container>::linear2(const typename grid<type,container>::pointer xout)
{
    assert(xout->get_dimension() == image_reference->get_dimension());

    int w = xout->get_size()[0];
    int h = xout->get_size()[1];

    // setup output image
    auto image_out = image<type,container>::new_pointer(w, h);
    image_out->set_sod_parameters(xout->get_spacing(), xout->get_origin(), xout->get_direction());
    image_out->assign(_default_);

    if (w < 1 || h < 1) return image_out;

    typename image<type,container>::pointer xo = (xout->ptr()[0]);   // raw pointer to
    typename image<type,container>::pointer yo = (xout->ptr()[1]);   // input grid coordinates

    // convert xout to x_reference (world to image coordinates)
    // std::cout << "convert xout to xref" << std::endl;
    std::vector<double> sod = image_reference->get_sod_inverse();
    auto xro = image<type,container>::new_pointer(xo->get_size());
    auto yro = image<type,container>::new_pointer(yo->get_size());
    container::affine_sod_2d(xo->get_data(), yo->get_data(), xro->get_data(), yro->get_data(), sod);

    // apply the interpolation algorithm (or kernel)
    // std::cout << "linear call to container" << std::endl;
    container::linear2(xro->get_data(),yro->get_data(),image_reference->get_data(),
                image_out->get_data(), image_reference->get_size(), xout->get_size() );
    // image_out->set_data(imgo);

    return image_out;
};

template <typename type, typename container>
typename image<type,container>::pointer ilinear<type,container>::linear3(const typename grid<type,container>::pointer xout)
{
    assert(xout->get_dimension() == image_reference->get_dimension());

    int w = xout->get_size()[0];
    int h = xout->get_size()[1];
    int l = xout->get_size()[2];

    // setup output image
    auto image_out = image<type,container>::new_pointer(w, h, l);
    image_out->set_sod_parameters(xout->get_spacing(), xout->get_origin(), xout->get_direction());
    image_out->assign(_default_);

    if (w < 1 || h < 1 || l < 1) return image_out;

    typename image<type,container>::pointer xo = (xout->ptr()[0]);   // raw pointer to
    typename image<type,container>::pointer yo = (xout->ptr()[1]);   // input grid coordinates
    typename image<type,container>::pointer zo = (xout->ptr()[2]);

    // convert xout to x_reference (world to image coordinates)
    // std::cout << "convert xout to xref" << std::endl;
    std::vector<double> sod = image_reference->get_sod_inverse();
    auto xro = image<type,container>::new_pointer(xo->get_size());
    auto yro = image<type,container>::new_pointer(yo->get_size());
    auto zro = image<type,container>::new_pointer(zo->get_size());
    container::affine_sod_3d(xo->get_data(), yo->get_data(), zo->get_data(),
                             xro->get_data(), yro->get_data(), zro->get_data(), sod);
    // xro->print_data();

    // apply the interpolation algorithm (or kernel)
    // std::cout << "linear call to container" << std::endl;
    container::linear3(xro->get_data(),yro->get_data(), zro->get_data(), image_reference->get_data(),
                image_out->get_data(), image_reference->get_size(), xout->get_size() );
    // image_out->set_data(imgo);
    // std::cout << "done" << std::endl;
    return image_out;
};

}; //end namespace

#endif