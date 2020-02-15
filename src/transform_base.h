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
public:
    //Type definitions
    using pointer = std::shared_ptr<transform_base<pixel_type>>;
    using vector = std::vector<transform_base::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    typename image<pixel_type>::pointer parameters;
    typename image<pixel_type>::pointer inverse_parameters;

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
    transform_base(int d, typename image<pixel_type>::pointer params);
    transform_base(const transform_base<pixel_type> & input);

    virtual void copy(const transform_base<pixel_type> & input);
    virtual void duplicate(const transform_base & input);

    // ===========================================
    // Get Functions
    // ===========================================
    typename image<pixel_type>::pointer get_parameters() const;
    typename image<pixel_type>::pointer get_inverse_parameters() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_parameters(typename image<pixel_type>::pointer params);

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
transform_base<pixel_type>::transform_base(int d, typename image<pixel_type>::pointer params)
{
    this->class_name = "transform_base";
    init(d);
    parameters = params;
    inverse_();
};

template <typename pixel_type>
void transform_base<pixel_type>::init(int d)
{
    typename image<pixel_type>::pointer param( std::make_shared<image<float>>(d) );
    typename image<pixel_type>::pointer inv( std::make_shared<image<float>>(d) );
    parameters = param;
    inverse_parameters = inv;
    
    object<pixel_type>::init(d);
};

template <typename pixel_type>
void transform_base<pixel_type>::copy(const transform_base<pixel_type> & input)
{
    object<pixel_type>::copy_properties(input);
    (*(this->parameters)).copy(*input.get_parameters());
    // (*(this->inverse_parameters)).copy(*input.get_inverse_parameters());
};


// Full copy
template <typename pixel_type>
void transform_base<pixel_type>::duplicate(const transform_base<pixel_type> & input)
{
    object<pixel_type>::copy_properties(input);

    this->parameters = input.get_parameters();
    this->inverse_parameters = input.get_inverse_parameters();
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
typename image<pixel_type>::pointer transform_base<pixel_type>::get_parameters() const
{
    // image<pixel_type> im;
    // std::cout << "a\n";
    // // parameters.print_data();
    // // im.copy(this->parameters);
    // std::cout << "a\n";
    // return im;
    return parameters;
};

template <typename pixel_type>
typename image<pixel_type>::pointer transform_base<pixel_type>::get_inverse_parameters() const
{
    return inverse_parameters;
};

// ===========================================
// Set Functions
// ===========================================
template <typename pixel_type>
void transform_base<pixel_type>::set_parameters(typename image<pixel_type>::pointer params)
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
    typename image<pixel_type>::pointer params = get_inverse_parameters();
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