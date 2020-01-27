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
#include "image.h"
#include "grid.h"

// Class transform_base
template <typename pixel_type>
class transform_base: public object<pixel_type>
{
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    image<pixel_type> parameters;
    image<pixel_type> inverse_parameters;

    // ===========================================
    // Functions
    // ===========================================
    virtual void init(int d);
    virtual void copy(const transform_base<pixel_type> & input);

    // Internal function (compute inverse)
    virtual void inverse_(); // parameters that do nothing when transform is applied

public:
    // ===========================================
    // Create Functions
    // ===========================================
    transform_base();
    transform_base(int d);
    transform_base(int d, image<pixel_type> & params);
    transform_base(const transform_base<pixel_type> & input);

    // ===========================================
    // Get Functions
    // ===========================================
    image<pixel_type> get_parameters() const;
    image<pixel_type> get_inverse_parameters() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_parameters(image<pixel_type> & params);

    // ===========================================
    // Print Functions
    // ===========================================
    // void print(std::string msg = "");
    std::string info(std::string msg = "");

    // template <typename pixel_type_>
    // friend std::ostream & operator << (std::ostream & os, transform_base<pixel_type_> & input);
    
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
    parameters = image<pixel_type>(2);
    init(2);
};

template <typename pixel_type>
transform_base<pixel_type>::transform_base(int d)
{
    this->class_name = "transform_base";
    parameters = image<pixel_type>(d);
    init(d);
};

template <typename pixel_type>
transform_base<pixel_type>::transform_base(int d, image<pixel_type> & params)
{
    this->class_name = "transform_base";
    parameters = params;
    init(d);
};

template <typename pixel_type>
void transform_base<pixel_type>::init(int d)
{
    inverse_();
    object<pixel_type>::init(d);
};

template <typename pixel_type>
void transform_base<pixel_type>::copy(const transform_base<pixel_type> & input)
{
    object<pixel_type>::copy_properties(input);
    image<pixel_type> params = input.get_parameters();
    this->parameters = params;
};

// ===========================================
// Get Functions
// ===========================================
template <typename pixel_type>
image<pixel_type> transform_base<pixel_type>::get_parameters() const
{
    return parameters;
};

template <typename pixel_type>
image<pixel_type> transform_base<pixel_type>::get_inverse_parameters() const
{
    return inverse_parameters;
};

// ===========================================
// Set Functions
// ===========================================
template <typename pixel_type>
void transform_base<pixel_type>::set_parameters(image<pixel_type> & params)
{
    parameters = params;
    inverse_();
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

    ss << object<pixel_type>::info(title);
    
    ss << "Parameters: \t\t";
    ss << parameters.info_data("");
    // ss << "]" << std::endl;
    // ss << std::endl;

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
    ;
};

template <typename pixel_type>
transform_base<pixel_type> transform_base<pixel_type>::inverse()
{
    image<pixel_type> params = get_inverse_parameters();
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
    return transform(point);
};

template <typename pixel_type>
grid<pixel_type> transform_base<pixel_type>::operator *(grid<pixel_type> & input)
{
    return transform(input);
};

#endif