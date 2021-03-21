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
#include "image_utils.h"
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
    // Type definitions
    using vector_image     = std::vector<typename image<type,container>::pointer>;
    using ptr_vector_image = std::shared_ptr<vector_image>;

    ptr_vector_image parameters;
    ptr_vector_image inverse_parameters;

    std::vector<double> sigma;
    std::vector<int> kernel;

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

    std::vector<double> get_sigma() const;
    std::vector<int> get_kernel() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_parameters(typename image<type,container>::pointer params, int n = 0);
    void set_parameters_vector(ptr_vector_image params);
    void set_sigma(std::vector<double> s);
    void set_kernel(std::vector<int> k);
    virtual void change_size(std::vector<int> sz);

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
    virtual typename transform<type,container>::pointer inverse(); // parameters that do nothing when transform is applied

    // ===========================================
    // Functions
    // ===========================================
    virtual std::vector<type> apply(std::vector<type> & point);
    // virtual grid<type,container> apply(const grid<type,container> & input);
    virtual typename grid<type,container>::pointer apply(typename grid<type,container>::pointer input);

    std::vector<type> operator () (std::vector<type> & point);
    // grid<type,container> operator () (const grid<type,container> & input);
    typename grid<type,container>::pointer operator () (typename grid<type,container>::pointer input);

    virtual void fluid();
    virtual void elastic();

    virtual void read(std::string file_name);
    virtual void write(std::string file_name);
};

template<typename type>
using transform_cpu = transform<type,vector_cpu<type>>;

// template<typename type>
// using transform_gpu = transform<type,vector_opencl<type>>;

#ifdef IMART_WITH_OPENCL
template<typename type>
using transform_opencl = transform<type,vector_opencl<type>>;
#endif

#ifdef IMART_WITH_CUDA
template<typename type>
using transform_cuda = transform<type,vector_cuda<type>>;
#endif

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
    // std::cout << "init allocate" << std::endl;
    parameters = std::make_shared< std::vector< typename image<type,container>::pointer >>(d);
    inverse_parameters = std::make_shared< std::vector< typename image<type,container>::pointer >>(d);
    for(int i = 0; i < d; i++)
    {
        (*parameters)[i] = image<type,container>::new_pointer();        // parameters are 2d
        (*inverse_parameters)[i] = image<type,container>::new_pointer();// parameters are 2d
    };
    // std::cout << "end allocate" << std::endl;
};

// Full copy
template <typename type, typename container>
void transform<type,container>::clone_(const transform<type,container> & input)
{
    space_object::mimic_(input);
    int n = input.get_parameters_vector()->size();
    allocate(n);
    for(int i = 0; i < n; i++)
    {
        (*parameters)[i]->clone_(*(input.get_parameters(i)));
        (*inverse_parameters)[i]->clone_(*(input.get_inverse_parameters(i)));
    };
    sigma = input.get_sigma();
    kernel = input.get_kernel();
};

template <typename type, typename container>
void transform<type,container>::copy_(const transform<type,container> & input)
{
    space_object::mimic_(input);
    int n = input.get_parameters_vector()->size();
    allocate(n);
    for(int i = 0; i < n; i++)
    {
        // std::cout << "copy" << std::endl;
        (*parameters)[i]->copy_(*(input.get_parameters(i)));
        (*inverse_parameters)[i]->copy_(*(input.get_inverse_parameters(i)));
    };
    sigma = input.get_sigma();
    kernel = input.get_kernel();
};

template <typename type, typename container>
void transform<type,container>::mimic_(const transform<type,container> & input)
{
    space_object::mimic_(input);
    int n = input.get_parameters_vector()->size();
    allocate(n);
    for(int i = 0; i < n; i++)
    {
        (*parameters)[i]->mimic_(*(input.get_parameters(i)));
        (*parameters)[i]->zeros();
        (*inverse_parameters)[i]->mimic_(*(input.get_inverse_parameters(i)));
    };
    sigma = input.get_sigma();
    kernel = input.get_kernel();
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

template <typename type, typename container>
std::vector<double> transform<type,container>::get_sigma() const
{
    return sigma;
};

template <typename type, typename container>
std::vector<int> transform<type,container>::get_kernel() const
{
    return kernel;
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

template <typename type, typename container>
void transform<type,container>::set_sigma(std::vector<double> s)
{
    sigma = s;
};

template <typename type, typename container>
void transform<type,container>::set_kernel(std::vector<int> k)
{
    kernel = k;
};

template <typename type, typename container>
void transform<type,container>::change_size(std::vector<int> sz)
{
    return;
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
    ss << "Size: \t\t\t[ ";
    for(int i = 0; i < this->size.size(); i++) { ss << this->size[i] << " "; };
    ss << "]" << std::endl;
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
typename transform<type,container>::pointer transform<type,container>::inverse()
{
    inverse_(); // update, compute in case of parameter update
    typename transform<type,container>::ptr_vector_image params = get_inverse_parameters_vector();
    // auto tr = transform<type,container>::new_pointer(this->dim, params);
    auto tr = this->mimic();
    tr->set_parameters_vector(params);
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

template <typename type, typename container>
void transform<type,container>::fluid()
{
    ;
};

template <typename type, typename container>
void transform<type,container>::elastic()
{
    ;
};

template <typename type, typename container>
void transform<type,container>::read(std::string file_name)
{
    ;
};

template <typename type, typename container>
void transform<type,container>::write(std::string file_name)
{
    ;
};

}; //end namespace

#endif