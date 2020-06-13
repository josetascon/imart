/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __VECTOR_GPU_H__
#define __vector_gpu_H__

// std libs
#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <cassert>      // assert

// gpu libs
#include <CL/cl.hpp>

// def
#ifndef VIENNACL_WITH_OPENCL
  #define VIENNACL_WITH_OPENCL
#endif

// viennacl headers
#include <viennacl/ocl/backend.hpp>
#include <viennacl/scalar.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/linalg/vector_operations.hpp>

// local libs
#include "inherit.h"
#include "object.h"
#include "kernels.h"

namespace imart
{

// Class object
template <typename type>
class vector_gpu: public inherit<vector_gpu<type>, object>, viennacl::vector<type>
{
public:
    //Type definitions
    using self    = vector_gpu;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using object::class_name;
    using object::get_name;
    // using viennacl::vector<type>::operator =;

    using viennacl::vector<type>::size;
    using viennacl::vector<type>::begin;
    using viennacl::vector<type>::end;
    using viennacl::vector<type>::operator[];
    using viennacl::vector<type>::handle;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    // int lower_limit;
    // int upper_limit;
    
public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    vector_gpu(): viennacl::vector<type>() { class_name = "vector_gpu"; };          // constructor empty
    vector_gpu(int s): viennacl::vector<type>(s) { class_name = "vector_gpu"; };    // constructor
    vector_gpu(int s, type value): viennacl::vector<type>(s) { class_name = "vector_gpu"; viennacl::linalg::vector_assign(*this, type(value)); };
    vector_gpu(const vector_gpu & input);       // constructor clone
    ~vector_gpu();                              // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const vector_gpu & input);  // copy everything
    virtual void copy_ (const vector_gpu & input);  // share data
    virtual void mimic_(const vector_gpu & input);  // copy meta data

    // ===========================================
    // Get Functions
    // ===========================================

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);
    std::string info_data(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    void assert_size(const vector_gpu<type> & input);

    // ===========================================
    // Initialize Functions
    // ===========================================
    // TODO
    // void zeros();
    // void ones();
    // void fill(type value);
    // void random(float min=0.0, float max=1.0);
    
    // ===========================================
    // Overloading operators
    // ===========================================
    // Vector to vector
    pointer operator = (const vector_gpu<type>::pointer input);
    pointer operator + (const vector_gpu<type> & input);
    pointer operator - (const vector_gpu<type> & input);
    pointer operator * (const vector_gpu<type> & input);
    pointer operator / (const vector_gpu<type> & input);
    pointer operator ^ (const vector_gpu<type> & input);

    // Scalar to vector
    pointer operator + (type scalar);
    pointer operator - (type scalar);
    pointer operator * (type scalar);
    pointer operator / (type scalar);
    pointer operator ^ (type scalar);

    // Friend classes to support scalar to vector left hand side
    template<typename type_>
    friend typename vector_gpu<type_>::pointer operator + (type_ scalar, vector_gpu<type_> & input);
    template<typename type_>
    friend typename vector_gpu<type_>::pointer operator - (type_ scalar, const vector_gpu<type_> & input);
    template<typename type_>
    friend typename vector_gpu<type_>::pointer operator * (type_ scalar, vector_gpu<type_> & input);
    template<typename type_>
    friend typename vector_gpu<type_>::pointer operator / (type_ scalar, const vector_gpu<type_> & input);

    // ===========================================
    // Reduction functions
    // ===========================================
    type min();
    type max();
    type sum();
    // type prod();     // may produce overflow error
    type dot(const vector_gpu<type> & input);

    // ===========================================
    // Functions
    // ===========================================
    pointer normalize(type min = 0.0, type max = 1.0);

    template<typename type_cast>
    typename vector_gpu<type_cast>::pointer cast();
};


// ===========================================
//          Functions of Class vector_gpu
// ===========================================

// ===========================================
// Constructors
// ===========================================
template <typename type>
vector_gpu<type>::vector_gpu(const vector_gpu<type> & input)
{
    clone_(input);          // call the virtual
};

//! Destructor
template <typename type>
vector_gpu<type>::~vector_gpu()
{
    ;
};

// ===========================================
// Create Functions
// ===========================================
template <typename type>
void vector_gpu<type>::clone_(const vector_gpu<type> & input)
{
    // std::cout << "clone_";
    mimic_(input);
    viennacl::copy(input.begin(),input.end(),this->begin());
};

template <typename type>
void vector_gpu<type>::copy_(const vector_gpu<type> & input)
{
    // std::cout << "copy_";
    mimic_(input);
};

template <typename type>
void vector_gpu<type>::mimic_(const vector_gpu<type> & input)
{
    // std::cout << "mimic_";
    int size = input.size();
    // this->resize(size, false); // 2nd parameter preserve; too slow!!
    // *this = vector_gpu<type>(size);  // equal; also slow!!
    this->~vector_gpu<type>();          // destruct
    new(this) vector_gpu<type>(size);   // reconstruct
};

// ===========================================
// Get Functions
// ===========================================

// ===========================================
// Print Functions
// ===========================================
template <typename type>
std::string vector_gpu<type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Vector GPU Information";
    if (msg != "") { title = msg; };

    // Summary of the object information
    ss << object::info(title);
    ss << "Size: \t\t\t" << this->size() << std::endl;
    // ss << "Capacity: \t\t" << this->capacity() << std::endl;
    return ss.str();
};

template <typename type>
std::string vector_gpu<type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; };
    
    int size = this->size();
    std::vector<type> tmp(size);    // create std::vector to copy and print, same as viennacl
    viennacl::copy(this->begin(), this->end(), tmp.begin());
    
    type * p = tmp.data();
    ss << "[ ";
    for(int k=0; k<tmp.size(); k++) { ss << p[k] << " "; };
    ss << "]" << std::endl;
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type>
void vector_gpu<type>::assert_size(const vector_gpu<type> & input)
{
    assert(this->size() == input.size());
};

// ===========================================
// Overloading operators
// ===========================================
template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator = (const vector_gpu<type>::pointer input)
{
    return input;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator + (const vector_gpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    // *output = viennacl::operator + ( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "+");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator - (const vector_gpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    // *output = viennacl::operator - ( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "-");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator * (const vector_gpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    // *output = viennacl::linalg::element_prod( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "*");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator / (const vector_gpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    // *output = viennacl::linalg::element_div( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "/");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator ^ (const vector_gpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    // *output = viennacl::linalg::element_pow( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "pow", true);
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

// Scalar right hand side
template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator + (type scalar)
{
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);

    std::string str_kernel = kernel_scalar( string_type<type>(), "+");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator - (type scalar)
{
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "-");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator * (type scalar)
{
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "*");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator / (type scalar)
{
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "/");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::operator ^ (type scalar)
{
    int size = this->size();
    vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "pow", true);  // function = true
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

// Scalar left hand side
template <typename type>
typename vector_gpu<type>::pointer operator + (type scalar, vector_gpu<type> & input)
{
    return input + scalar;
};

template <typename type>
typename vector_gpu<type>::pointer operator - (type scalar, const vector_gpu<type> & input)
{
    int size = input.size();
    typename vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "-", false, true);  // function = true, reverse = true
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(input, *output, scalar));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer operator * (type scalar, vector_gpu<type> & input)
{
    return input * scalar;
};

template <typename type>
typename vector_gpu<type>::pointer operator / (type scalar, const vector_gpu<type> & input)
{
    int size = input.size();
    typename vector_gpu<type>::pointer output = vector_gpu<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "/", false, true);  // function = true, reverse = true
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(input, *output, scalar));
    return output;
};

// ===========================================
// Reduction Functions
// ===========================================
template <typename type>
type vector_gpu<type>::min()
{
    type x;
    viennacl::scalar<type> xm(0);
    viennacl::linalg::min_impl(*this, xm);
    x = xm;
    return x;
};

template <typename type>
type vector_gpu<type>::max()
{
    type x;
    viennacl::scalar<type> xm(0);
    viennacl::linalg::max_impl(*this, xm);
    x = xm;
    return x;
};

template <typename type>
type vector_gpu<type>::sum()
{
    type x;
    viennacl::scalar<type> xm(0);
    viennacl::linalg::sum_impl(*this, xm);
    x = xm;
    return x;
};

// Vectorial dot product. Verify the same number of elements, then product and reduce
template <typename type>
type vector_gpu<type>::dot(const vector_gpu<type> & input)
{
    assert_size(input);

    type x;
    viennacl::scalar<type> xm(0);
    viennacl::linalg::inner_prod_impl(*this, input, xm);
    x = xm;
    return x;
};

// ===========================================
// Functions
// ===========================================
template <typename type> template <typename type_cast>
typename vector_gpu<type_cast>::pointer vector_gpu<type>::cast()
{
    int size = this->size();
    typename vector_gpu<type_cast>::pointer output = vector_gpu<type_cast>::new_pointer(size);
    // viennacl::linalg::convert(*output,*this);
    std::string str_kernel = kernel_cast( string_type<type>(), string_type<type_cast>());
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_cast");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_cast");
    viennacl::ocl::enqueue(kkk(*this, *output));
    return output;
};

template <typename type>
typename vector_gpu<type>::pointer vector_gpu<type>::normalize(type min, type max)
{
    type minv = this->min();
    type maxv = this->max();

    vector_gpu<type>::pointer output = *this - minv;
    output = *(*output*((max - min)/(maxv - minv))) + min;
    return output->clone(); // memory error if viennacl is out of context, clone for keeping alive!!!
};

}; //end namespace

#endif