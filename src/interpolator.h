/*
* @Author: jose
* @Date:   2019-11-19 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-19 13:55:00
*/

#ifndef __INTERPOLATOR_H__
#define __INTERPOLATOR_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector

// local libs
#include "process_object.h"
#include "image.h"
#include "grid.h"

namespace imart
{

// Class interpolator
template <typename type, typename container=vector_cpu<type>>
class interpolator: public inherit<interpolator<type,container>, process_object>
{
public:
    //Type definitions
    using self    = interpolator;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using process_object::init;
    using process_object::info;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    typename image<type,container>::pointer image_reference;
    typename image<type,container>::pointer image_output;
    typename grid<type,container>::pointer x_output;
    // typename grid<type,container>::pointer x_reference;
    // borders<type> region;
    type _default_;

    // ===========================================
    // Functions
    // ===========================================
    void init(int d);
    
public:
    // ===========================================
    // Constructor Function
    // ===========================================
    interpolator();
    interpolator(int d);
    interpolator(typename image<type,container>::pointer imgref);
    // interpolator(typename image<type,container>::pointer imgref, typename grid<type,container>::pointer xref);

    // ===========================================
    // Create Functions
    // ===========================================
    void clone_(const interpolator<type,container> & input);
    void copy_(const interpolator<type,container> & input);
    void mimic_(const interpolator<type,container> & input);

    // ===========================================
    // Get Functions
    // ===========================================
    type get_default() const;
    typename image<type,container>::pointer get_reference() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_default(type value);
    void set_reference(typename image<type,container>::pointer imgref);

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    virtual typename image<type,container>::pointer apply(const typename grid<type,container>::pointer xout);
    typename image<type,container>::pointer operator () (const typename grid<type,container>::pointer xout);
};


// ===========================================
//      Functions of Class grid
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
template <typename type, typename container>
interpolator<type,container>::interpolator()
{
    this->class_name = "interpolator";
    init(2);
    // image_reference = imgref;
    // x_reference = xref;
};

template <typename type, typename container>
interpolator<type,container>::interpolator(int d)
{
    this->class_name = "interpolator";
    init(d);
    // image_reference = imgref;
    // x_reference = xref;
};

// template <typename type, typename container>
// interpolator<type,container>::interpolator(typename image<type,container>::pointer imgref, typename grid<type,container>::pointer xref)
template <typename type, typename container>
interpolator<type,container>::interpolator(typename image<type,container>::pointer imgref)
{
    // assert(imgref->get_dimension() == xref->get_dimension());
    this->class_name = "interpolator";
    int d = imgref->get_dimension();
    init(d);
    image_reference = imgref;
    // x_reference = xref;
};

template <typename type, typename container>
void interpolator<type,container>::init(int d)
{
    _default_ = 0;                  // fill values that are not interpolator
    this->set_total_inputs(2);      //process_object::init
    this->set_total_outputs(1);     //process_object::init
    // this->set_total_inputs(3);      //process_object::init using x_reference
    // std::cout << "inputs = " << this->num_inputs << std::endl; 
    // std::cout << "outputs = " << this->num_outputs << std::endl; 
    
    image_reference = image<type,container>::new_pointer(d);
    x_output = grid<type,container>::new_pointer(d);
    image_output = image<type,container>::new_pointer(d);
    // x_reference = grid<type,container>::new_pointer(d);

    // this->setup_input(image_reference, x_reference, x_output);
    this->setup_input(image_reference, x_output);
    this->setup_output(image_output);
};

// ===========================================
// Create Functions
// ===========================================
template <typename type, typename container>
void interpolator<type,container>::clone_(const interpolator<type,container> & input)
{
    process_object::clone_(input);
    image_reference->clone_(*(input.get_reference()));
};

template <typename type, typename container>
void interpolator<type,container>::copy_(const interpolator<type,container> & input)
{
    process_object::copy_(input);
    image_reference->copy_(*(input.get_reference()));
};

template <typename type, typename container>
void interpolator<type,container>::mimic_(const interpolator<type,container> & input)
{
    process_object::copy_(input); // the mimic do no copy the input or output vectors
    image_reference->copy_(*(input.get_reference())); // copy the image (eficient)
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
type interpolator<type,container>::get_default() const
{
    return _default_;
};

template <typename type, typename container>
typename image<type,container>::pointer interpolator<type,container>::get_reference() const
{
    return image_reference;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void interpolator<type,container>::set_default(type value)
{
    _default_ = value;
};

template <typename type, typename container>
void interpolator<type,container>::set_reference(typename image<type,container>::pointer imgref)
{
    image_reference = imgref;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string interpolator<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Interpolator Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    ss << process_object::info("");
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
typename image<type,container>::pointer interpolator<type,container>::apply (const typename grid<type,container>::pointer xout)
{
    return image<type,container>::new_pointer();
}


template <typename type, typename container>
typename image<type,container>::pointer interpolator<type,container>::operator () (const typename grid<type,container>::pointer xout)
{
    return apply(xout);
}

}; //end namespace


#endif