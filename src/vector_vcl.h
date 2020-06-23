/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __VECTOR_VCL_H__
#define __VECTOR_VCL_H__

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
class vector_vcl: public inherit<vector_vcl<type>, object>, viennacl::vector<type>
{
public:
    //Type definitions
    using self    = vector_vcl;
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
    vector_vcl(): viennacl::vector<type>() { class_name = "vector_vcl"; };          // constructor empty
    vector_vcl(int s): viennacl::vector<type>(s) { class_name = "vector_vcl"; };    // constructor
    vector_vcl(int s, type value): viennacl::vector<type>(s) { class_name = "vector_vcl"; viennacl::linalg::vector_assign(*this, type(value)); };
    vector_vcl(std::initializer_list<type> list);
    vector_vcl(const vector_vcl & input);       // constructor clone
    ~vector_vcl();                              // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const vector_vcl & input);  // copy everything
    virtual void copy_ (const vector_vcl & input);  // share data
    virtual void mimic_(const vector_vcl & input);  // copy meta data

    // ===========================================
    // Get Functions
    // ===========================================
    std::vector<type> std_vector();

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    void assert_size(const vector_vcl<type> & input);

    // ===========================================
    // Initialize Functions
    // ===========================================
    void zeros();
    void ones();
    void assign(type value);
    void random(float min=0.0, float max=1.0);
    
    // ===========================================
    // Overloading operators
    // ===========================================
    // Vector to vector
    pointer operator = (const vector_vcl<type>::pointer input);
    pointer operator + (const vector_vcl<type> & input);
    pointer operator - (const vector_vcl<type> & input);
    pointer operator * (const vector_vcl<type> & input);
    pointer operator / (const vector_vcl<type> & input);
    pointer operator ^ (const vector_vcl<type> & input);

    // Scalar to vector
    pointer operator + (type scalar);
    pointer operator - (type scalar);
    pointer operator * (type scalar);
    pointer operator / (type scalar);
    pointer operator ^ (type scalar);

    // Friend classes to support scalar to vector left hand side
    template<typename type_>
    friend typename vector_vcl<type_>::pointer operator + (type_ scalar, vector_vcl<type_> & input);
    template<typename type_>
    friend typename vector_vcl<type_>::pointer operator - (type_ scalar, const vector_vcl<type_> & input);
    template<typename type_>
    friend typename vector_vcl<type_>::pointer operator * (type_ scalar, vector_vcl<type_> & input);
    template<typename type_>
    friend typename vector_vcl<type_>::pointer operator / (type_ scalar, const vector_vcl<type_> & input);

    // ===========================================
    // Reduction functions
    // ===========================================
    type min();
    type max();
    type sum();
    // type prod();     // may produce overflow error
    type dot(const vector_vcl<type> & input);

    // ===========================================
    // Functions
    // ===========================================
    pointer normalize(type min = 0.0, type max = 1.0);

    template<typename type_cast>
    typename vector_vcl<type_cast>::pointer cast();

    static std::vector<typename vector_vcl<type>::pointer> grid_2d(int w, int h, std::vector<double> & sod);
    static std::vector<typename vector_vcl<type>::pointer> grid_3d(int w, int h, int l, std::vector<double> & sod);

};


// ===========================================
//          Functions of Class vector_vcl
// ===========================================

// ===========================================
// Constructors
// ===========================================
template <typename type>
vector_vcl<type>::vector_vcl(std::initializer_list<type> list)
{
    class_name = "vector_vcl";
    int size = list.size();
    this->~vector_vcl<type>();          // destruct
    new(this) vector_vcl<type>(size);   // reconstruct
    viennacl::copy(list.begin(),list.end(),this->begin()); //check if works!
};

template <typename type>
vector_vcl<type>::vector_vcl(const vector_vcl<type> & input)
{
    class_name = "vector_vcl";
    clone_(input);          // call the virtual
};

//! Destructor
template <typename type>
vector_vcl<type>::~vector_vcl()
{
    ;
};

// ===========================================
// Create Functions
// ===========================================
template <typename type>
void vector_vcl<type>::clone_(const vector_vcl<type> & input)
{
    // std::cout << "clone_";
    mimic_(input);
    viennacl::copy(input.begin(),input.end(),this->begin());
};

template <typename type>
void vector_vcl<type>::copy_(const vector_vcl<type> & input)
{
    // std::cout << "copy_";
    mimic_(input);
};

template <typename type>
void vector_vcl<type>::mimic_(const vector_vcl<type> & input)
{
    // std::cout << "mimic_";
    int size = input.size();
    // this->resize(size, false); // 2nd parameter preserve; too slow!!
    // *this = vector_vcl<type>(size);  // equal; also slow!!
    this->~vector_vcl<type>();          // destruct
    new(this) vector_vcl<type>(size);   // reconstruct
};

// ===========================================
// Get Functions
// ===========================================
template <typename type>
std::vector<type> vector_vcl<type>::std_vector()
{
    std::vector<type> tmp(this->size());
    viennacl::copy(this->begin(), this->end(), tmp.begin());
    return tmp;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type>
std::string vector_vcl<type>::info(std::string msg)
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
std::string vector_vcl<type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; };

    std::vector<type> tmp = std_vector();
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
void vector_vcl<type>::assert_size(const vector_vcl<type> & input)
{
    assert(this->size() == input.size());
};

// ===========================================
// Initialization Functions
// ===========================================
template <typename type>
void vector_vcl<type>::zeros()
{
    viennacl::linalg::vector_assign(*this, type(0));
};

template <typename type>
void vector_vcl<type>::ones()
{
    viennacl::linalg::vector_assign(*this, type(1));
};

// Set all pixel to a fixed value
template <typename type>
void vector_vcl<type>::assign(type value)
{
    viennacl::linalg::vector_assign(*this, type(value));
};

template <typename type>
void vector_vcl<type>::random(float min, float max)
{
    std::string str_kernel = kernel_random( string_type<type>(), min, max);
    // std::cout << str_kernel << std::endl;
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_random");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_random");
    viennacl::ocl::enqueue(kkk(*this));

    // std::string str_kernel2 = kernel_random( string_type<type>(), min, max, false);
    // std::cout << str_kernel2 << std::endl;
    // viennacl::ocl::program & prog2 = viennacl::ocl::current_context().add_program(str_kernel2, "kernel_random_next");
    // viennacl::ocl::kernel & kkk2 = prog2.get_kernel("kernel_random_next");
    // viennacl::ocl::enqueue(kkk2(*this));
};

// ===========================================
// Overloading operators
// ===========================================
template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator = (const vector_vcl<type>::pointer input)
{
    return input;
};

template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator + (const vector_vcl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);
    // *output = viennacl::operator + ( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "+");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator - (const vector_vcl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);
    // *output = viennacl::operator - ( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "-");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator * (const vector_vcl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);
    // *output = viennacl::linalg::element_prod( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "*");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator / (const vector_vcl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);
    // *output = viennacl::linalg::element_div( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "/");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator ^ (const vector_vcl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);
    // *output = viennacl::linalg::element_pow( *this, input );
    std::string str_kernel = kernel_vector( string_type<type>(), "pow", true);
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_vector");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_vector");
    viennacl::ocl::enqueue(kkk(*this, input, *output));
    return output;
};

// Scalar right hand side
template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator + (type scalar)
{
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);

    std::string str_kernel = kernel_scalar( string_type<type>(), "+");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator - (type scalar)
{
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "-");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator * (type scalar)
{
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "*");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator / (type scalar)
{
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "/");
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::operator ^ (type scalar)
{
    int size = this->size();
    auto output = vector_vcl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "pow", true);  // function = true
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(*this, *output, scalar));
    return output;
};

// Scalar left hand side
template <typename type>
typename vector_vcl<type>::pointer operator + (type scalar, vector_vcl<type> & input)
{
    return input + scalar;
};

template <typename type>
typename vector_vcl<type>::pointer operator - (type scalar, const vector_vcl<type> & input)
{
    int size = input.size();
    auto output = vector_vcl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "-", false, true);  // function = true, reverse = true
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_scalar");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_scalar");
    viennacl::ocl::enqueue(kkk(input, *output, scalar));
    return output;
};

template <typename type>
typename vector_vcl<type>::pointer operator * (type scalar, vector_vcl<type> & input)
{
    return input * scalar;
};

template <typename type>
typename vector_vcl<type>::pointer operator / (type scalar, const vector_vcl<type> & input)
{
    int size = input.size();
    auto output = vector_vcl<type>::new_pointer(size);
    
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
type vector_vcl<type>::min()
{
    type x;
    viennacl::scalar<type> xm(0);
    viennacl::linalg::min_impl(*this, xm);
    x = xm;
    return x;
};

template <typename type>
type vector_vcl<type>::max()
{
    type x;
    viennacl::scalar<type> xm(0);
    viennacl::linalg::max_impl(*this, xm);
    x = xm;
    return x;
};

template <typename type>
type vector_vcl<type>::sum()
{
    type x;
    viennacl::scalar<type> xm(0);
    viennacl::linalg::sum_impl(*this, xm);
    x = xm;
    return x;
};

// Vectorial dot product. Verify the same number of elements, then product and reduce
template <typename type>
type vector_vcl<type>::dot(const vector_vcl<type> & input)
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
template <typename type>
typename vector_vcl<type>::pointer vector_vcl<type>::normalize(type min, type max)
{
    type minv = this->min();
    type maxv = this->max();

    vector_vcl<type>::pointer output = *this - minv;
    output = *(*output*((max - min)/(maxv - minv))) + min;
    return output->clone(); // memory error if viennacl is out of context, clone for keeping alive!!!
};

template <typename type> template <typename type_cast>
typename vector_vcl<type_cast>::pointer vector_vcl<type>::cast()
{
    int size = this->size();
    auto output = vector_vcl<type_cast>::new_pointer(size);
    // viennacl::linalg::convert(*output,*this);
    std::string str_kernel = kernel_cast( string_type<type>(), string_type<type_cast>());
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_cast");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_cast");
    viennacl::ocl::enqueue(kkk(*this, *output));
    return output;
};

// template <typename type>
// std::vector<typename vector_vcl<type>::pointer> grid_2d(int w, int h)
template <typename type>
std::vector<typename vector_vcl<type>::pointer> vector_vcl<type>::grid_2d(int w, int h, std::vector<double> & sod)
{
    int size = w*h;
    auto x = vector_vcl<type>::new_pointer(size);
    auto y = vector_vcl<type>::new_pointer(size);
    
    auto p = vector_vcl<double>::new_pointer(sod.size());
    viennacl::copy(sod.begin(),sod.end(),p->begin());

    std::string str_kernel = kernel_grid_2d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_grid_2d");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_grid_2d");
    viennacl::ocl::enqueue(kkk(*x, *y, *p, w, h));
    
    std::vector<typename vector_vcl<type>::pointer> xy(2);
    xy[0] = x;
    xy[1] = y;
    return xy;
};

template <typename type>
std::vector<typename vector_vcl<type>::pointer> vector_vcl<type>::grid_3d(int w, int h, int l, std::vector<double> & sod)
{
    int size = w*h*l;
    auto x = vector_vcl<type>::new_pointer(size);
    auto y = vector_vcl<type>::new_pointer(size);
    auto z = vector_vcl<type>::new_pointer(size);

    auto p = vector_vcl<double>::new_pointer(sod.size());
    viennacl::copy(sod.begin(),sod.end(),p->begin());

    std::string str_kernel = kernel_grid_3d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    viennacl::ocl::program & prog = viennacl::ocl::current_context().add_program(str_kernel, "kernel_grid_3d");
    viennacl::ocl::kernel & kkk = prog.get_kernel("kernel_grid_3d");
    viennacl::ocl::enqueue(kkk(*x, *y, *z, *p, w, h, l));
    
    std::vector<typename vector_vcl<type>::pointer> xyz(3);
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    return xyz;
};

}; //end namespace

#endif