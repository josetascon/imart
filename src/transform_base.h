/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-15 13:55:00
*/

#ifndef __TRANSFORM_BASE_H__
#define __TRANSFORM_BASE_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <cassert>      // assert

// local libs 
#include "object.h"
#include "image_base.h"
#include "grid.h"

namespace imart
{

// Class transform_base
template <typename pixel_type>
class transform_base: public object
{
public:
    //Type definitions
    using self    = transform_base;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    typename image_base<pixel_type>::pointer parameters;
    typename image_base<pixel_type>::pointer inverse_parameters;

    // ===========================================
    // Functions
    // ===========================================
    virtual void init(int d);

    // Internal function (compute inverse)
    virtual void inverse_(); // parameters that do nothing when transform is applied

public:
    // ===========================================
    // Create Functions
    // ===========================================
    transform_base();
    transform_base(int d);
    transform_base(int d, typename image_base<pixel_type>::pointer params);
    transform_base(const transform_base<pixel_type> & input);

    template<typename... ARGS>
    static pointer new_pointer(const ARGS&... args);

    void copy(const transform_base<pixel_type> & input);
    void duplicate(const transform_base & input);
    void imitate(const transform_base & input);

    // ===========================================
    // Get Functions
    // ===========================================
    typename image_base<pixel_type>::pointer get_parameters() const;
    typename image_base<pixel_type>::pointer get_inverse_parameters() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_parameters(typename image_base<pixel_type>::pointer params);

    // ===========================================
    // Print Functions
    // ===========================================
    // void print(std::string msg = "");
    std::string info(std::string msg);
    std::string info_data(std::string msg);

    // template <typename pixel_type_>
    // friend std::ostream & operator << (std::ostream & os, transform_base<pixel_type_> & input);
    
    // ===========================================
    // Overloading Functions
    // ===========================================
    // Equal
    transform_base<pixel_type> & operator = (const transform_base<pixel_type> & input);

    // ===========================================
    // Initialization Functions
    // ===========================================
    virtual void identity(); // parameters that do nothing when transform is applied
    virtual transform_base<pixel_type> inverse(); // parameters that do nothing when transform is applied

    // ===========================================
    // Functions
    // ===========================================
    virtual std::vector<pixel_type> transform(std::vector<pixel_type> & point);
    virtual grid<pixel_type> transform(grid<pixel_type> & input);

    std::vector<pixel_type> operator * (std::vector<pixel_type> & point);
    grid<pixel_type> operator * (grid<pixel_type> & input);

    transform_base<pixel_type> operator + (const transform_base<pixel_type> & input);
    transform_base<pixel_type> operator - (const transform_base<pixel_type> & input);
    transform_base<pixel_type> operator * (pixel_type scalar);
    transform_base<pixel_type> operator * (const image_base<pixel_type> & input);
};


// ===========================================
//          Functions of Class transform_base
// ===========================================

// ===========================================
// Create Functions
// ===========================================
template <typename pixel_type>
transform_base<pixel_type>::transform_base()
{
    this->class_name = "transform_base";
    init(2);
};

template <typename pixel_type>
transform_base<pixel_type>::transform_base(int d)
{
    this->class_name = "transform_base";
    init(d);
};

template <typename pixel_type>
transform_base<pixel_type>::transform_base(int d, typename image_base<pixel_type>::pointer params)
{
    this->class_name = "transform_base";
    init(d);
    parameters = params;
    inverse_();
};

template <typename pixel_type>
transform_base<pixel_type>::transform_base(const transform_base<pixel_type> & input)
{
    this->class_name = "transform_base";
    copy(input);
};

template <typename pixel_type>
void transform_base<pixel_type>::init(int d)
{
    // typename image_base<pixel_type>::pointer param( std::make_shared<image_base<pixel_type>>(d) );
    // typename image_base<pixel_type>::pointer inv( std::make_shared<image_base<pixel_type>>(d) );
    // parameters = param;
    // inverse_parameters = inv;
    parameters = image_base<pixel_type>::new_pointer(d);
    inverse_parameters = image_base<pixel_type>::new_pointer(d);
    
    object::init(d);
};

template <typename pixel_type>
template <typename ... ARGS>
typename transform_base<pixel_type>::pointer transform_base<pixel_type>::new_pointer(const ARGS&... args)
{
    return std::make_shared< transform_base<pixel_type> >(args...); // not working for inherited classes
};

// Full copy
template <typename pixel_type>
void transform_base<pixel_type>::copy(const transform_base<pixel_type> & input)
{
    object::copy_properties(input);
    // (*(parameters)).copy(*input.get_parameters());
    // (*(inverse_parameters)).copy(*input.get_inverse_parameters());

    parameters->copy(*input.get_parameters());
    inverse_parameters->copy(*input.get_inverse_parameters());
};

template <typename pixel_type>
void transform_base<pixel_type>::duplicate(const transform_base<pixel_type> & input)
{
    object::copy_properties(input);

    parameters = input.get_parameters();
    inverse_parameters = input.get_inverse_parameters();
};

template <typename pixel_type>
void transform_base<pixel_type>::imitate(const transform_base<pixel_type> & input)
{
    object::copy_properties(input);
    parameters->imitate(*input.get_parameters());
    inverse_parameters->imitate(*input.get_inverse_parameters());
};

// Equal
template <typename pixel_type>
transform_base<pixel_type> & transform_base<pixel_type>::operator = (const transform_base<pixel_type> & input)
{
    // delete &data;
    duplicate(input);
    return *this;
};

// ===========================================
// Get Functions
// ===========================================
template <typename pixel_type>
typename image_base<pixel_type>::pointer transform_base<pixel_type>::get_parameters() const
{
    return parameters;
};

template <typename pixel_type>
typename image_base<pixel_type>::pointer transform_base<pixel_type>::get_inverse_parameters() const
{
    return inverse_parameters;
};

// ===========================================
// Set Functions
// ===========================================
template <typename pixel_type>
void transform_base<pixel_type>::set_parameters(typename image_base<pixel_type>::pointer params)
{
    parameters = params;
    // inverse_();
};



// ===========================================
// Print Functions
// ===========================================
// template <typename pixel_type>
// void transform_base<pixel_type>::print(std::string msg)
// {
//     std::cout << transform_base::info(msg);
// };

template <typename pixel_type>
std::string transform_base<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Transform Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    
    ss << "Number parameters: \t";
    ss << parameters->get_total_elements() << std::endl;
    // ss << "]" << std::endl;
    // ss << std::endl;

    return ss.str();
};

template <typename pixel_type>
std::string transform_base<pixel_type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; }
    else { ss << "Transform parameters:\n"; };
    ss << parameters->info_data("");

    return ss.str();
};


// template <typename pixel_type>
// std::ostream & operator << (std::ostream & os, transform_base<pixel_type> & input)
// {
//     os << input.info("");
//     return os;
// };


// ===========================================
// Initialization Functions
// ===========================================
template <typename pixel_type>
void transform_base<pixel_type>::identity()
{
    ;
};

template <typename pixel_type>
void transform_base<pixel_type>::inverse_()
{
    (*inverse_parameters).copy(*parameters);
};

template <typename pixel_type>
transform_base<pixel_type> transform_base<pixel_type>::inverse()
{
    typename image_base<pixel_type>::pointer params = get_inverse_parameters();
    transform_base<pixel_type> tr(this->dim, params);
    return tr;
};

// ===========================================
// Functions
// ===========================================
template <typename pixel_type>
std::vector<pixel_type> transform_base<pixel_type>::transform(std::vector<pixel_type> & point)
{
    return point;
};

template <typename pixel_type>
grid<pixel_type> transform_base<pixel_type>::transform(grid<pixel_type> & input)
{
    return input;
};

template <typename pixel_type>
std::vector<pixel_type> transform_base<pixel_type>::operator *(std::vector<pixel_type> & point)
{
    return this->transform(point);
};

template <typename pixel_type>
grid<pixel_type> transform_base<pixel_type>::operator *(grid<pixel_type> & input)
{
    return this->transform(input);
};

// Transform to Transform
template <typename pixel_type>
transform_base<pixel_type> transform_base<pixel_type>::operator + (const transform_base<pixel_type> & input)
{
    auto pp = image_base<pixel_type>::new_pointer();
    *pp = *parameters + *(input.get_parameters());

    transform_base<pixel_type> result;
    result.imitate(*this);
    result.set_parameters( pp );
    return result;
};

template <typename pixel_type>
transform_base<pixel_type> transform_base<pixel_type>::operator - (const transform_base<pixel_type> & input)
{
    auto pp = image_base<pixel_type>::new_pointer();
    *pp = *parameters - *(input.get_parameters());

    transform_base<pixel_type> result;
    result.imitate(*this);
    result.set_parameters( pp );
    return result;
};

template <typename pixel_type>
transform_base<pixel_type> transform_base<pixel_type>::operator * (pixel_type scalar)
{
    auto pp = image_base<pixel_type>::new_pointer();
    *pp = scalar*(*parameters);

    transform_base<pixel_type> result;
    result.imitate(*this);
    result.set_parameters( pp );
    return result;
};

template <typename pixel_type>
transform_base<pixel_type> transform_base<pixel_type>::operator * (const image_base<pixel_type> & input)
{
    auto pp = image_base<pixel_type>::new_pointer();
    *pp = input*(*parameters);

    transform_base<pixel_type> result;
    result.imitate(*this);
    result.set_parameters( pp );
    return result;
};

}; //end namespace

#endif