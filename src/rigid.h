/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __RIGID_H__
#define __RIGID_H__

// local libs
#include "global.h"

namespace imart
{

// Class rigid
template <typename type, typename container=vector_cpu<type>>
class rigid: public inherit<rigid<type,container>, global<type,container>>
{
public:
    //Type definitions
    using self    = rigid;
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

    using inherit<rigid<type,container>, global<type,container>>::inherit;

protected:
    // ===========================================
    // Functions
    // ===========================================
    // void init(int d);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    rigid()      : inherit<rigid<type,container>, global<type,container>>()  { this->class_name = "rigid"; };
    rigid(int d) : inherit<rigid<type,container>, global<type,container>>(d) { this->class_name = "rigid"; };
    rigid(int d, typename image<type,container>::pinter params) 
                 : inherit<rigid<type,container>, global<type,container>>(d, params) { this->class_name = "rigid"; };
    
    rigid( type thetaz, type tx, type ty );
    rigid( type thetaz, type thetay, type thetax, type tx, type ty, type tz );
};

template<typename type>
using rigid_cpu = rigid<type,vector_cpu<type>>;

template<typename type>
using rigid_gpu = rigid<type,vector_ocl<type>>;

template<typename type>
using rigid_cuda = rigid<type,vector_cuda<type>>;


// ===========================================
//          Functions of Class rigid
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
template <typename type, typename container>
rigid<type,container>::rigid( type thetaz, type tx, type ty )
{
    this->class_name = "rigid";
    init(2);
    std::initializer_list<type> ll{ cos(thetaz), -sin(thetaz), sin(thetaz), cos(thetaz), tx, ty };
    this->set_parameters( image<type,container>::new_pointer(ll) );
    inverse_();
};

template <typename type, typename container>
rigid<type,container>::rigid( type thetaz, type thetay, type thetax, type tx, type ty, type tz )
{
    this->class_name = "rigid";
    init(3);
    std::initializer_list<type> ll{ cos(thetaz)*cos(thetay), cos(thetaz)*sin(thetay)*sin(thetax) - sin(thetaz)*cos(thetax),  sin(thetaz)*sin(thetax) - cos(thetaz)*sin(thetay)*cos(thetax),
                                    sin(thetaz)*cos(thetay), sin(thetaz)*sin(thetay)*sin(thetax) + cos(thetaz)*cos(thetax), -sin(thetaz)*sin(thetay)*cos(thetax) - cos(thetaz)*sin(thetax),
                                    -sin(thetay)           , cos(thetay)*sin(thetax)                                      ,  cos(thetay)*cos(thetax),    
                                    tx, ty, tz };
    this->set_parameters( image<type,container>::new_pointer(ll) );
    inverse_();
};


}; //end namespace

#endif