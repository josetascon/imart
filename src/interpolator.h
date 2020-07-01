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

template <typename type>
struct borders
{
    type xmin;
    type xmax;
    type ymin;
    type ymax;
    type zmin;
    type zmax;
};

// Class interpolator
template <typename type, typename container=vector_cpu<type>>
class interpolator: public inherit<interpolator<type,container>, process_object>
{
public:
    //Type definitions
    using self    = interpolator;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    using process_object::init;
    using process_object::info;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    typename image<type,container>::pointer image_reference;
    typename image<type,container>::pointer image_output;
    typename grid<type,container>::pointer x_reference;
    typename grid<type,container>::pointer x_output;
    borders<type> region;
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
    interpolator(typename image<type,container>::pointer imgref, typename grid<type,container>::pointer xref);

    // ===========================================
    // Get Functions
    // ===========================================
    type get_default() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_default(type value);

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
interpolator<type,container>::interpolator(typename image<type,container>::pointer imgref, typename grid<type,container>::pointer xref)
{
    assert(imgref->get_dimension() == xref->get_dimension());
    
    this->class_name = "interpolator";
    int d = imgref->get_dimension();
    init(d);
    image_reference = imgref;
    x_reference = xref;
};

template <typename type, typename container>
void interpolator<type,container>::init(int d)
{
    this->set_total_inputs(3);      //process_object::init
    this->set_total_outputs(1);     //process_object::init         
    // std::cout << "inputs = " << this->num_inputs << std::endl; 
    // std::cout << "outputs = " << this->num_outputs << std::endl; 
    _default_ = 0;                  // fill values that are not interpolatord
    
    image_reference = image<type,container>::new_pointer(d);
    x_reference = grid<type,container>::new_pointer(d);
    x_output = grid<type,container>::new_pointer(d);
    image_output = image<type,container>::new_pointer(d);

    this->setup_input(image_reference, x_reference, x_output);
    this->setup_output(image_output);
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
type interpolator<type,container>::get_default() const
{
    return _default_;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void interpolator<type,container>::set_default(type value)
{
    _default_ = value;
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