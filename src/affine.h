/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __AFFINE_H__
#define __AFFINE_H__

// local libs
#include "global.h"

namespace imart
{

// Class affine
template <typename type, typename container=vector_cpu<type>>
class affine: public inherit<affine<type,container>, global<type,container>>
{
public:
    //Type definitions
    using self    = affine;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using transform<type,container>::parameters;
    using transform<type,container>::inverse_parameters;
    using transform<type,container>::get_parameters;
    using transform<type,container>::get_inverse_parameters;
    using transform<type,container>::operator=;

    using global<type,container>::init;
    using global<type,container>::inverse_;

    using inherit<affine<type,container>, global<type,container>>::inherit;

protected:
    // ===========================================
    // Functions
    // ===========================================
    // void init(int d);
    
public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    affine()      : inherit<affine<type,container>, global<type,container>>()  { this->class_name = "affine"; };
    affine(int d) : inherit<affine<type,container>, global<type,container>>(d) { this->class_name = "affine"; };
    affine(int d, typename image<type,container>::pointer params)
                  : inherit<affine<type,container>, global<type,container>>(d, params) { this->class_name = "affine"; };

    affine( type a1, type a2, type a3, type a4, type tx, type ty );
    affine( type a1, type a2, type a3, type a4, type a5, type a6, type a7, type a8, type a9, type tx, type ty, type tz );

};

template<typename type>
using affine_cpu = affine<type,vector_cpu<type>>;

// template<typename type>
// using affine_gpu = affine<type,vector_opencl<type>>;

#ifdef IMART_WITH_OPENCL
template<typename type>
using affine_opencl = affine<type,vector_opencl<type>>;
#endif

#ifdef IMART_WITH_CUDA
template<typename type>
using affine_cuda = affine<type,vector_cuda<type>>;
#endif

// ===========================================
//          Functions of Class affine
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
template <typename type, typename container>
affine<type,container>::affine( type a1, type a2, type a3, type a4, type tx, type ty )
{
    this->class_name = "affine";
    init(2);
    std::initializer_list<type> ll{ a1, a2, a3, a4, tx, ty };
    this->set_parameters( image<type,container>::new_pointer(ll) );
    inverse_();
};

template <typename type, typename container>
affine<type,container>::affine( type a1, type a2, type a3, type a4, type a5, type a6, type a7, type a8, type a9, type tx, type ty, type tz )
{
    this->class_name = "affine";
    init(3);
    std::initializer_list<type> ll{ a1, a2, a3, a4, a5, a6, a7, a8, a9, tx, ty, tz };
    this->set_parameters( image<type,container>::new_pointer(ll) );
    inverse_();
};

}; //end namespace

#endif