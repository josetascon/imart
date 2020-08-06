/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-15 13:55:00
*/

#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <cassert>      // assert

// local libs
#include "space_object.h"
#include "image.h"
#include "grid.h"

namespace imart
{

// Class transform
template <typename type, typename container=vector_cpu<type>>
class transform: public inherit<transform<type,container>, space_object>
{
public:
    //Type definitions
    using self    = transform;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    typename image<type,container>::pointer parameters;
    typename image<type,container>::pointer inverse_parameters;

    // ===========================================
    // Functions
    // ===========================================
    virtual void init(int d);
    virtual void inverse_(); // Internal function (compute inverse). Parameters that do nothing when transform is applied

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    transform();
    transform(int d);
    transform(int d, typename image<type,container>::pointer params);
    transform(const transform<type,container> & input);

    // ===========================================
    // Create Functions
    // ===========================================
    void clone_(const transform & input);
    void copy_(const transform & input);
    void mimic_(const transform & input);

    // ===========================================
    // Get Functions
    // ===========================================
    typename image<type,container>::pointer get_parameters() const;
    typename image<type,container>::pointer get_inverse_parameters() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_parameters(typename image<type,container>::pointer params);

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);
    std::string info_data(std::string msg);

    // ===========================================
    // Overloading Functions
    // ===========================================
    // transform<type,container> & operator = (const transform<type,container> & input);
    // transform<type,container> operator + (const transform<type,container> & input);
    // transform<type,container> operator - (const transform<type,container> & input);
    // transform<type,container> operator * (const image<type,container> & input);
    // transform<type,container> operator * (type scalar);

    // ===========================================
    // Initialization Functions
    // ===========================================
    virtual void identity(); // parameters that do nothing when transform is applied
    virtual transform<type,container> inverse(); // parameters that do nothing when transform is applied

    // ===========================================
    // Functions
    // ===========================================
    virtual std::vector<type> apply(std::vector<type> & point);
    // virtual grid<type,container> apply(const grid<type,container> & input);
    virtual typename grid<type,container>::pointer apply(typename grid<type,container>::pointer input);

    std::vector<type> operator () (std::vector<type> & point);
    // grid<type,container> operator () (const grid<type,container> & input);
    typename grid<type,container>::pointer operator () (typename grid<type,container>::pointer input);
};


// ===========================================
//          Functions of Class transform
// ===========================================

// ===========================================
// Create Functions
// ===========================================
template <typename type, typename container>
transform<type,container>::transform()
{
    this->class_name = "transform";
    init(2);
};

template <typename type, typename container>
transform<type,container>::transform(int d)
{
    this->class_name = "transform";
    init(d);
};

template <typename type, typename container>
transform<type,container>::transform(int d, typename image<type,container>::pointer params)
{
    this->class_name = "transform";
    init(d);
    parameters = params;
    inverse_();
};

template <typename type, typename container>
transform<type,container>::transform(const transform<type,container> & input)
{
    this->class_name = "transform";
    clone_(input);
};

template <typename type, typename container>
void transform<type,container>::init(int d)
{
    // typename image<type,container>::pointer param( std::make_shared<image<type,container>>(d) );
    // typename image<type,container>::pointer inv( std::make_shared<image<type,container>>(d) );
    // parameters = param;
    // inverse_parameters = inv;

    parameters = image<type,container>::new_pointer();          // parameters are 2d
    inverse_parameters = image<type,container>::new_pointer();  // parameters are 2d
    space_object::init(d);
    // identity();
};

// Full copy
template <typename type, typename container>
void transform<type,container>::clone_(const transform<type,container> & input)
{
    space_object::mimic_(input);
    parameters->clone_(*(input.get_parameters()));
    inverse_parameters->clone_(*(input.get_inverse_parameters()));
};

template <typename type, typename container>
void transform<type,container>::copy_(const transform<type,container> & input)
{
    space_object::mimic_(input);
    parameters = input.get_parameters();
    inverse_parameters = input.get_inverse_parameters();
};

template <typename type, typename container>
void transform<type,container>::mimic_(const transform<type,container> & input)
{
    space_object::mimic_(input);
    parameters->mimic_(*(input.get_parameters()));
    inverse_parameters->mimic_(*(input.get_inverse_parameters()));
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
typename image<type,container>::pointer transform<type,container>::get_parameters() const
{
    return parameters;
};

template <typename type, typename container>
typename image<type,container>::pointer transform<type,container>::get_inverse_parameters() const
{
    return inverse_parameters;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void transform<type,container>::set_parameters(typename image<type,container>::pointer params)
{
    parameters = params;
    // inverse_();
};



// ===========================================
// Print Functions
// ===========================================
// template <typename type, typename container>
// void transform<type,container>::print(std::string msg)
// {
//     std::cout << transform::info(msg);
// };

template <typename type, typename container>
std::string transform<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Transform Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    ss << "Data type: \t\t" << parameters->get_type() << std::endl;
    ss << "Number parameters: \t" << parameters->get_total_elements() << std::endl;
    ss << "Parameters pointer: \t" << parameters << std::endl;
    ss << "Inverse pointer: \t" << inverse_parameters << std::endl;
    return ss.str();
};

template <typename type, typename container>
std::string transform<type,container>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; }
    else { ss << "Transform parameters:\n"; };
    ss << parameters->info_data("");

    return ss.str();
};


// template <typename type, typename container>
// std::ostream & operator << (std::ostream & os, transform<type,container> & input)
// {
//     os << input.info("");
//     return os;
// };


// ===========================================
// Initialization Functions
// ===========================================
template <typename type, typename container>
void transform<type,container>::identity()
{
    ;
};

template <typename type, typename container>
void transform<type,container>::inverse_()
{
    inverse_parameters->clone_(*parameters); // in base class just copy
};

template <typename type, typename container>
transform<type,container> transform<type,container>::inverse()
{
    inverse_(); // update, compute in case of parameter update
    typename image<type,container>::pointer params = get_inverse_parameters();
    transform<type,container> tr(this->dim, params);
    return tr;
};

// ===========================================
// Overloading Functions
// ===========================================
/*
// Equal
template <typename type, typename container>
transform<type,container> & transform<type,container>::operator = (const transform<type,container> & input)
{
    // delete &data;
    copy_(input);
    return *this;
};

// Transform to Transform
template <typename type, typename container>
transform<type,container> transform<type,container>::operator + (const transform<type,container> & input)
{
    auto pp = image<type,container>::new_pointer();
    *pp = *parameters + *(input.get_parameters());

    pointer output = this->mimic();
    output->set_parameters( pp );
    return *output;
    // transform<type,container> output;
    // output.mimic_(input);
    // output.set_parameters( pp );
    // return output;
};

template <typename type, typename container>
transform<type,container> transform<type,container>::operator - (const transform<type,container> & input)
{
    auto pp = image<type,container>::new_pointer();
    *pp = *parameters - *(input.get_parameters());

    pointer output = this->mimic();
    output->set_parameters( pp );
    return *output;
};

template <typename type, typename container>
transform<type,container> transform<type,container>::operator * (const image<type,container> & input)
{
    auto pp = image<type,container>::new_pointer();
    *pp = input*(*parameters);
    pointer output = this->mimic();
    output->set_parameters( pp );
    return *output;
};

// Scalar
template <typename type, typename container>
transform<type,container> transform<type,container>::operator * (type scalar)
{
    auto pp = image<type,container>::new_pointer();
    *pp = scalar*(*parameters);

    pointer output = this->mimic();
    output->set_parameters( pp );
    return *output;
};*/

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
std::vector<type> transform<type,container>::apply(std::vector<type> & point)
{
    return point;
};

// template <typename type, typename container>
// grid<type,container> transform<type,container>::apply(const grid<type,container> & input)
// {
//     return input;
// };

template <typename type, typename container>
typename grid<type,container>::pointer transform<type,container>::apply(typename grid<type,container>::pointer input)
{
    return input;
};

template <typename type, typename container>
std::vector<type> transform<type,container>::operator() (std::vector<type> & point)
{
    return apply(point);
};

// template <typename type, typename container>
// grid<type,container> transform<type,container>::operator() (const grid<type,container> & input)
// {
//     return apply(input);
// };

template <typename type, typename container>
typename grid<type,container>::pointer transform<type,container>::operator() (typename grid<type,container>::pointer input)
{
    return apply(input);
};

}; //end namespace

#endif