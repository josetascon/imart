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
#include "image_utils.h"

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
    // Type definitions
    using vector_image     = std::vector<typename image<type,container>::pointer>;
    using ptr_vector_image = std::shared_ptr<vector_image>;

    ptr_vector_image parameters;
    ptr_vector_image inverse_parameters;

    // ===========================================
    // Functions
    // ===========================================
    virtual void init(int d);
    virtual void allocate(int d);
    virtual void inverse_(); // Internal function (compute inverse). Parameters that do nothing when transform is applied

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    transform();
    transform(int d);
    transform(int d, typename image<type,container>::pointer params);
    transform(int d, typename transform<type,container>::ptr_vector_image params);
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
    typename image<type,container>::pointer get_parameters(int n = 0) const;
    typename image<type,container>::pointer get_inverse_parameters(int n = 0) const;

    ptr_vector_image get_parameters_vector() const;
    ptr_vector_image get_inverse_parameters_vector() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_parameters(typename image<type,container>::pointer params, int n = 0);
    void set_parameters_vector(ptr_vector_image params);

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);
    std::string info_data(std::string msg);

    // ===========================================
    // Overloading Functions
    // ===========================================
    transform<type,container> & operator = (const transform<type,container> & input);
    transform<type,container> operator + (const transform<type,container> & input);
    transform<type,container> operator - (const transform<type,container> & input);
    transform<type,container> operator * (const transform<type,container> & input);
    transform<type,container> operator * (const image<type,container> & input);
    transform<type,container> operator * (type scalar);

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
    (*parameters)[0] = params;
    inverse_();
};

template <typename type, typename container>
transform<type,container>::transform(int d, typename transform<type,container>::ptr_vector_image params)
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
    // parameters = image<type,container>::new_pointer();          // parameters are 2d
    // inverse_parameters = image<type,container>::new_pointer();  // parameters are 2d
    // space_object::init(d);
    // // identity();
    allocate(d);
    space_object::init(d);
};

template <typename type, typename container>
void transform<type,container>::allocate(int d)
{
    parameters = std::make_shared< std::vector< typename image<type,container>::pointer >>(d);
    inverse_parameters = std::make_shared< std::vector< typename image<type,container>::pointer >>(d);
    for(int i = 0; i < d; i++)
    {
        (*parameters)[i] = image<type,container>::new_pointer();        // parameters are 2d
        (*inverse_parameters)[i] = image<type,container>::new_pointer();// parameters are 2d
    };
};

// Full copy
template <typename type, typename container>
void transform<type,container>::clone_(const transform<type,container> & input)
{
    int n = input.get_parameters_vector()->size();
    allocate(n);
    space_object::mimic_(input);
    for(int i = 0; i < n; i++)
    {
        (*parameters)[i]->clone_(*(input.get_parameters(i)));
        (*inverse_parameters)[i]->clone_(*(input.get_inverse_parameters(i)));
    };
};

template <typename type, typename container>
void transform<type,container>::copy_(const transform<type,container> & input)
{
    int n = input.get_parameters_vector()->size();
    allocate(n);
    space_object::mimic_(input);
    for(int i = 0; i < n; i++)
    {
        // std::cout << "copy" << std::endl;
        (*parameters)[i]->copy_(*(input.get_parameters(i)));
        (*inverse_parameters)[i]->copy_(*(input.get_inverse_parameters(i)));
    };
};

template <typename type, typename container>
void transform<type,container>::mimic_(const transform<type,container> & input)
{
    int n = input.get_parameters_vector()->size();
    allocate(n);
    space_object::mimic_(input);
    for(int i = 0; i < n; i++)
    {
        (*parameters)[i]->mimic_(*(input.get_parameters(i)));
        (*inverse_parameters)[i]->mimic_(*(input.get_inverse_parameters(i)));
    };
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
typename image<type,container>::pointer transform<type,container>::get_parameters(int n) const
{
    return (*parameters)[n];
};

template <typename type, typename container>
typename image<type,container>::pointer transform<type,container>::get_inverse_parameters(int n) const
{
    return (*inverse_parameters)[n];
};

template <typename type, typename container>
typename transform<type,container>::ptr_vector_image transform<type,container>::get_parameters_vector() const
{
    return parameters;
};

template <typename type, typename container>
typename transform<type,container>::ptr_vector_image transform<type,container>::get_inverse_parameters_vector() const
{
    return inverse_parameters;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void transform<type,container>::set_parameters(typename image<type,container>::pointer params, int n)
{
    (*parameters)[n] = params;
    // inverse_();
};

template <typename type, typename container>
void transform<type,container>::set_parameters_vector(typename transform<type,container>::ptr_vector_image params)
{
    parameters = params;
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
    ss << "Data type: \t\t" << (*parameters)[0]->get_type() << std::endl;
    ss << "Group of parameters: \t" << parameters->size() << std::endl;
    ss << "Number parameters: \t" << (*parameters)[0]->get_total_elements() << std::endl;
    ss << "Parameters pointer: \t" << parameters;// << std::endl;
    for(int i = 0; i < parameters->size(); i++)
        ss << ", " << (*parameters)[i];
    ss << std::endl;
    ss << "Inverse pointer: \t" << inverse_parameters;// << std::endl;
    for(int i = 0; i < inverse_parameters->size(); i++)
        ss << ", " << (*inverse_parameters)[i];
    ss << std::endl;
    return ss.str();
};

template <typename type, typename container>
std::string transform<type,container>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; }
    else { ss << "Transform parameters:\n"; };
    
    for(int i = 0; i < parameters->size(); i++)
        ss << (*parameters)[i]->info_data("");

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
    int n = (*parameters).size();
    for(int i = 0; i < n; i++)
        (*inverse_parameters)[i]->clone_( *((*parameters)[i]) ); // in base class just copy
};

template <typename type, typename container>
transform<type,container> transform<type,container>::inverse()
{
    inverse_(); // update, compute in case of parameter update
    typename transform<type,container>::ptr_vector_image params = get_inverse_parameters_vector();
    transform<type,container> tr(this->dim, params);
    return tr;
};

// ===========================================
// Overloading Functions
// ===========================================

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
    // std::cout << "add transform" << std::endl;
    auto pp = std::make_shared< std::vector< typename image<type,container>::pointer >>(get_parameters_vector()->size());
    *pp = *parameters + *(input.get_parameters_vector());

    // pointer output = this->mimic();
    // // output->set_parameters( pp );
    // return *output;
    transform<type,container> output;
    output.mimic_(*this);
    output.set_parameters_vector( pp );
    return output;
};

template <typename type, typename container>
transform<type,container> transform<type,container>::operator - (const transform<type,container> & input)
{
    auto pp = std::make_shared< std::vector< typename image<type,container>::pointer >>(get_parameters_vector()->size());
    *pp = *parameters - *(input.get_parameters_vector());

    transform<type,container> output;
    output.mimic_(*this);
    output.set_parameters_vector( pp );
    return output;
};

template <typename type, typename container>
transform<type,container> transform<type,container>::operator * (const transform<type,container> & input)
{
    auto pp = std::make_shared< std::vector< typename image<type,container>::pointer >>(get_parameters_vector()->size());
    *pp = (*parameters)*(*input.get_parameters_vector());

    transform<type,container> output;
    output.mimic_(*this);
    output.set_parameters_vector( pp );
    return output;
};

// Scalar
template <typename type, typename container>
transform<type,container> transform<type,container>::operator * (type scalar)
{
    auto pp = std::make_shared< std::vector< typename image<type,container>::pointer >>(get_parameters_vector()->size());
    *pp = (*parameters)*scalar;

    transform<type,container> output;
    output.mimic_(*this);
    output.set_parameters_vector( pp );
    return output;
};

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