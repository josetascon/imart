/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __VECTOR_OCL_H__
#define __VECTOR_OCL_H__

// std libs
#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <random>       // std::random
#include <cassert>      // assert

// opencl libs
#include <CL/cl.hpp>
#include <clFFT.h>

// local libs
#include "object.h"
#include "kernels.h"
#include "opencl_object.h"
#include "utils/type.h"
#include "utils/timer.h"

namespace imart
{

opencl_object cl_manager; // manager of vector ocl;

// imart function to check opencl device
void imart_opencl_device_name()
{
    cl_manager.print_device_name();
};

// Class object
template <typename type>
class vector_ocl: public inherit<vector_ocl<type>, object>
{
public:
    //Type definitions
    using self    = vector_ocl;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using object::class_name;
    using object::get_name;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int _size_;
    cl_int err;
    std::shared_ptr<cl::Buffer> buffer;

    // ===========================================
    // Constructor Functions
    // ===========================================
    void init(int s);
    void allocate(int s);

    // ===========================================
    // Reduction functions
    // ===========================================
    pointer sum_loop(pointer input, cl_int workGroupSize, cl_int numWorkGroups);
    pointer minmax_loop(pointer input, cl_int workGroupSize, cl_int numWorkGroups);
    
public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    vector_ocl();                               // constructor empty
    vector_ocl(int s);                          // constructor size
    vector_ocl(int s, type value);              // constructor assign
    vector_ocl(std::initializer_list<type> list);
    vector_ocl(const vector_ocl & input);       // constructor clone
    ~vector_ocl();                              // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const vector_ocl & input);  // copy everything
    virtual void copy_ (const vector_ocl & input);  // share data
    virtual void mimic_(const vector_ocl & input);  // copy meta data

    // ===========================================
    // Get Functions
    // ===========================================
    int size() const;
    std::shared_ptr<cl::Buffer> get_buffer() const;
    std::vector<type> std_vector();

    // ===========================================
    // Memory Functions
    // ===========================================
    void read_ram(type * p, int size, int offset = 0);
    void write_ram(type * p, int size, int offset = 0);

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    void assert_size(const vector_ocl<type> & input);

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
    // Read
    type operator [](int e);

    // Vector to vector
    vector_ocl<type> & operator = (const vector_ocl & input);
    pointer operator = (const vector_ocl<type>::pointer input);
    pointer operator + (const vector_ocl<type> & input);
    pointer operator - (const vector_ocl<type> & input);
    pointer operator * (const vector_ocl<type> & input);
    pointer operator / (const vector_ocl<type> & input);
    pointer operator ^ (const vector_ocl<type> & input);

    // // Scalar to vector
    pointer operator + (type scalar);
    pointer operator - (type scalar);
    pointer operator * (type scalar);
    pointer operator / (type scalar);
    pointer operator ^ (type scalar);

    // Friend classes to support scalar to vector left hand side
    template<typename type_>
    friend typename vector_ocl<type_>::pointer operator + (type_ scalar, vector_ocl<type_> & input);
    template<typename type_>
    friend typename vector_ocl<type_>::pointer operator - (type_ scalar, const vector_ocl<type_> & input);
    template<typename type_>
    friend typename vector_ocl<type_>::pointer operator * (type_ scalar, vector_ocl<type_> & input);
    template<typename type_>
    friend typename vector_ocl<type_>::pointer operator / (type_ scalar, const vector_ocl<type_> & input);
    template<typename type_>
    friend typename vector_ocl<type_>::pointer operator ^ (type_ scalar, const vector_ocl<type_> & input);

    // Conditions 
    pointer operator == (const vector_ocl<type> & input);
    pointer operator > (const vector_ocl<type> & input);
    pointer operator < (const vector_ocl<type> & input);
    pointer operator >= (const vector_ocl<type> & input);
    pointer operator <= (const vector_ocl<type> & input);

    pointer operator == (type scalar);
    pointer operator > (type scalar);
    pointer operator < (type scalar); // std::shared_ptr<vector_ocl<type>>
    pointer operator >= (type scalar);
    pointer operator <= (type scalar);

    void replace(const vector_ocl<type> & idxs, const vector_ocl<type> & input);
    void replace(const vector_ocl<type> & idxs, type value);

    // ===========================================
    // Reduction functions
    // ===========================================
    type min();
    type max();
    type sum();
    type dot(const vector_ocl<type> & input);
    // type prod();     // may produce overflow error

    // // ===========================================
    // // Functions
    // // ===========================================
    pointer normalize(type min = 0.0, type max = 1.0);

    template<typename type_cast>
    typename vector_ocl<type_cast>::pointer cast();

    static void pad(typename vector_ocl<type>::pointer input, typename vector_ocl<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post);
    static void unpad(typename vector_ocl<type>::pointer input, typename vector_ocl<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post);

    static void affine_2d(typename vector_ocl<type>::pointer xi, typename vector_ocl<type>::pointer yi,
                          typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo,
                          typename vector_ocl<type>::pointer p);

    static void affine_3d(typename vector_ocl<type>::pointer xi, typename vector_ocl<type>::pointer yi, typename vector_ocl<type>::pointer zi,
                          typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo, typename vector_ocl<type>::pointer zo,
                          typename vector_ocl<type>::pointer p);

    static void affine_sod_2d(typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo,
                              typename vector_ocl<type>::pointer xr, typename vector_ocl<type>::pointer yr, 
                              std::vector<double> & sod);
    static void affine_sod_3d(typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo, typename vector_ocl<type>::pointer zo,
                              typename vector_ocl<type>::pointer xr, typename vector_ocl<type>::pointer yr, typename vector_ocl<type>::pointer zr,
                              std::vector<double> & sod);

    static void dfield_2d(typename vector_ocl<type>::pointer xi, typename vector_ocl<type>::pointer yi,
                          typename vector_ocl<type>::pointer x, typename vector_ocl<type>::pointer y,
                          typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo);

    static void dfield_3d(typename vector_ocl<type>::pointer xi, typename vector_ocl<type>::pointer yi, typename vector_ocl<type>::pointer zi,
                          typename vector_ocl<type>::pointer x, typename vector_ocl<type>::pointer y, typename vector_ocl<type>::pointer z,
                          typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo, typename vector_ocl<type>::pointer zo);

    static std::vector<typename vector_ocl<type>::pointer> grid2(int w, int h, std::vector<double> & sod);
    static std::vector<typename vector_ocl<type>::pointer> grid3(int w, int h, int l, std::vector<double> & sod);

    static void nearest2( typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo,
                        typename vector_ocl<type>::pointer imgr, typename vector_ocl<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);
    static void nearest3( typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo, typename vector_ocl<type>::pointer zo,
                        typename vector_ocl<type>::pointer imgr, typename vector_ocl<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);

    static void linear2( typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo,
                        typename vector_ocl<type>::pointer imgr, typename vector_ocl<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);
    static void linear3( typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo, typename vector_ocl<type>::pointer zo,
                        typename vector_ocl<type>::pointer imgr, typename vector_ocl<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);

    static void cubic2( typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo,
                        typename vector_ocl<type>::pointer imgr, typename vector_ocl<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);
    static void cubic3( typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo, typename vector_ocl<type>::pointer zo,
                        typename vector_ocl<type>::pointer imgr, typename vector_ocl<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);

    static void fft(std::vector<pointer> & input, std::vector<pointer> & output, std::vector<int> size, bool forward);

    static void gradientx( typename vector_ocl<type>::pointer imgr,
                           typename vector_ocl<type>::pointer imgo,
                           std::vector<int> ref_size);

    static void gradienty( typename vector_ocl<type>::pointer imgr,
                           typename vector_ocl<type>::pointer imgo,
                           std::vector<int> ref_size);

    static void gradientz( typename vector_ocl<type>::pointer imgr,
                           typename vector_ocl<type>::pointer imgo,
                           std::vector<int> ref_size);

    static void convolution( typename vector_ocl<type>::pointer imgr,
                             typename vector_ocl<type>::pointer kernel,
                             typename vector_ocl<type>::pointer imgo,
                             std::vector<int> ref_size, int kwidth);
};


// ===========================================
//          Functions of Class vector_ocl
// ===========================================

// ===========================================
// Constructors
// ===========================================
template <typename type>
vector_ocl<type>::vector_ocl()
{
    class_name = "vector_ocl";
    init(0);
};

template <typename type>
vector_ocl<type>::vector_ocl(int s)
{
    class_name = "vector_ocl";
    init(s);
};

template <typename type>
vector_ocl<type>::vector_ocl(int s, type value)
{
    class_name = "vector_ocl";
    init(s);
    assign(value);
};

template <typename type>
vector_ocl<type>::vector_ocl(std::initializer_list<type> list)
{
    class_name = "vector_ocl";
    int s = list.size();
    init(s);
    if (s > 0)
        cl_manager.get_queue().enqueueWriteBuffer(*buffer, CL_TRUE, 0, sizeof(type)*(s), list.begin());
};

template <typename type>
vector_ocl<type>::vector_ocl(const vector_ocl<type> & input)
{
    class_name = "vector_ocl";
    clone_(input);          // call the virtual
};

//! Destructor
template <typename type>
vector_ocl<type>::~vector_ocl()
{
    ;
};

template <typename type>
void vector_ocl<type>::init(int s)
{
    _size_ = s;
    if(_size_>0) allocate(_size_);
    else buffer = nullptr;
};

template <typename type>
void vector_ocl<type>::allocate(int s)
{
    err = 0;
    buffer = std::make_shared<cl::Buffer>(cl_manager.get_context(), CL_MEM_READ_WRITE, sizeof(type)*s, nullptr, &err);
    assert(err == 0);
    // zeros(); // init for gpu
};

// ===========================================
// Create Functions
// ===========================================
template <typename type>
void vector_ocl<type>::clone_(const vector_ocl<type> & input)
{
    // std::cout << "clone_";
    mimic_(input);
    // clone buffer;
    std::string str_kernel = kernel_copy(string_type<type>());
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_copy");
    cl_manager.arguments(*(input.get_buffer()), *buffer);
    cl_manager.execute(_size_);
};

template <typename type>
void vector_ocl<type>::copy_(const vector_ocl<type> & input)
{
    // std::cout << "copy_";
    _size_ = input.size();
    buffer.reset();
    buffer = input.get_buffer();
};

template <typename type>
void vector_ocl<type>::mimic_(const vector_ocl<type> & input)
{
    // std::cout << "mimic_";
    _size_ = input.size();
    buffer.reset();
    allocate(_size_);
};

// ===========================================
// Get Functions
// ===========================================
template <typename type>
int vector_ocl<type>::size() const
{
    return _size_;
};

template <typename type>
std::shared_ptr<cl::Buffer> vector_ocl<type>::get_buffer() const
{
    return buffer;
};

template <typename type>
std::vector<type> vector_ocl<type>::std_vector()
{
    std::vector<type> tmp(_size_);
    if (_size_ > 0)
        cl_manager.get_queue().enqueueReadBuffer(*buffer, CL_TRUE, 0, sizeof(type)*tmp.size(), tmp.data());
    return tmp;
};

// ===========================================
// Memory Functions
// ===========================================
template <typename type>
void vector_ocl<type>::read_ram(type * p, int s, int offset)
{
    if(s > 0)
        cl_manager.get_queue().enqueueWriteBuffer(*buffer, CL_TRUE, sizeof(type)*offset, sizeof(type)*(s), p);
};

template <typename type>
void vector_ocl<type>::write_ram(type * p, int s, int offset)
{
    if(s > 0)
        cl_manager.get_queue().enqueueReadBuffer(*buffer, CL_TRUE, sizeof(type)*offset, sizeof(type)*s, p);
};


// ===========================================
// Print Functions
// ===========================================
template <typename type>
std::string vector_ocl<type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Vector OpenCL Information";
    if (msg != "") { title = msg; };

    // Summary of the object information
    ss << object::info(title);
    ss << "Size: \t\t\t" << _size_ << std::endl;
    ss << "Memory: \t\t" << buffer << std::endl;
    // ss << "Capacity: \t\t" << this->capacity() << std::endl;
    return ss.str();
};

template <typename type>
std::string vector_ocl<type>::info_data(std::string msg)
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
void vector_ocl<type>::assert_size(const vector_ocl<type> & input)
{
    assert(this->size() == input.size());
};

// ===========================================
// Initialization Functions
// ===========================================
template <typename type>
void vector_ocl<type>::zeros()
{
    assign(type(0));
};

template <typename type>
void vector_ocl<type>::ones()
{
    assign(type(1));
};

// Set all pixel to a fixed value
template <typename type>
void vector_ocl<type>::assign(type value)
{
    type v = value;
    std::string str_kernel = kernel_assign(string_type<type>());
    cl_manager.program(str_kernel, "kernel_assign");
    cl_manager.arguments(*buffer, v);
    cl_manager.execute(_size_);
};

template <typename type>
void vector_ocl<type>::random(float min, float max)
{
    // std::string str_kernel = kernel_random( string_type<type>(), min, max);
    // // // std::cout << str_kernel << std::endl;
    // cl_manager.program(str_kernel, "kernel_random");
    // cl_manager.arguments(*buffer);
    // cl_manager.execute(_size_);

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    // std::default_random_engine gen(rd()); //Standard random generator()
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> uniform(min, max);
    int size = this->size();
    std::vector<type> vec(size);
    type * p = vec.data();
    
    for(int k=0; k<size; k++)
    {
        p[k] = (type)uniform(gen); // casting to pixel_type
    };

    read_ram(vec.data(),size);
};

// ===========================================
// Overloading operators
// ===========================================
template <typename type>
type vector_ocl<type>::operator [](int e)
{
    std::vector<type> tmp(1);
    if(e >= 0 && e < _size_ && _size_ > 0)
        // std::cout << "cl read" << std::endl;
        cl_manager.get_queue().enqueueReadBuffer(*buffer, CL_TRUE, sizeof(type)*(e), sizeof(type)*(1), tmp.data());
    return tmp[0];
};

template <typename type>
vector_ocl<type> & vector_ocl<type>::operator = (const vector_ocl & input)
{
    copy_(input);
    return *this;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator = (const vector_ocl<type>::pointer input)
{
    return input;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator + (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_vector( string_type<type>(), "+");
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator - (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_vector( string_type<type>(), "-");
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator * (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_vector( string_type<type>(), "*");
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator / (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_vector( string_type<type>(), "/");
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator ^ (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_vector( string_type<type>(), "pow", true);
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

// Scalar right hand side
template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator + (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_scalar( string_type<type>(), "+");
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator - (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "-");
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator * (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "*");
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator / (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "/");
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator ^ (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "pow", true);  // function = true
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

// Scalar left hand side
template <typename type>
typename vector_ocl<type>::pointer operator + (type scalar, vector_ocl<type> & input)
{
    return input + scalar;
};

template <typename type>
typename vector_ocl<type>::pointer operator - (type scalar, const vector_ocl<type> & input)
{
    int size = input.size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "-", false, true);  // function = false, reverse = true
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*(input.get_buffer()), *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer operator * (type scalar, vector_ocl<type> & input)
{
    return input * scalar;
};

template <typename type>
typename vector_ocl<type>::pointer operator / (type scalar, const vector_ocl<type> & input)
{
    int size = input.size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "/", false, true);  // function = false, reverse = true
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*(input.get_buffer()), *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer operator ^ (type scalar, const vector_ocl<type> & input)
{
    int size = input.size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "pow", true, true);  // function = false, reverse = true
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*(input.get_buffer()), *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator == (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_vector( string_type<type>(), "==");
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator > (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_vector( string_type<type>(), ">");
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator < (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_vector( string_type<type>(), "<");
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator >= (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_vector( string_type<type>(), ">=");
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator <= (const vector_ocl<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);

    std::string str_kernel = kernel_vector( string_type<type>(), "<=");
    cl_manager.program(str_kernel, "kernel_vector");
    cl_manager.arguments(*buffer, *(input.get_buffer()), *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator == (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "==", false);  // function = true
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator > (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), ">", false);  // function = true
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator < (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "<", false);  // function = true
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator >= (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), ">=", false);  // function = true
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::operator <= (type scalar)
{
    int size = this->size();
    auto output = vector_ocl<type>::new_pointer(size);
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "<=", false);  // function = true
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*buffer, *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

template <typename type>
void vector_ocl<type>::replace(const vector_ocl<type> & idxs, const vector_ocl<type> & input)
{
    assert_size(idxs);
    assert_size(input);
    int size = this->size();

    std::string str_kernel = kernel_replace( string_type<type>(), false);
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_replace");
    cl_manager.arguments(*(idxs.get_buffer()), *buffer, *(input.get_buffer()));
    cl_manager.execute(size);
    return;
};

template <typename type>
void vector_ocl<type>::replace(const vector_ocl<type> & idxs, type value)
{
    assert_size(idxs);
    int size = this->size();

    std::string str_kernel = kernel_replace( string_type<type>(), true);
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_replace");
    cl_manager.arguments( *(idxs.get_buffer()), *buffer, value );
    cl_manager.execute(size);
    return;
};

// ===========================================
// Reduction Functions
// ===========================================
template <typename type>
type vector_ocl<type>::min()
{
    // Pointer to work
    pointer input = this->copy();
    pointer output = input;

    // Working groups
    auto workGroupSize = cl_manager.get_work_group_size();
    auto numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    // std::cout << "\nwork group size =" << workGroupSize << std::endl;
    // std::cout << "num groups =" << numWorkGroups << std::endl;

    // Setup the kernel
    std::string str_kernel = kernel_minmax( string_type<type>(), false );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_min");
    
    while(numWorkGroups > workGroupSize)
    {
        output = minmax_loop(input, workGroupSize, numWorkGroups);
        input = output;
        numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    };

    // Finish the reduction work in CPU
    std::vector<type> out = output->std_vector();
    type min_ = out[0];
    for(int i = 0; i < out.size(); i++) min_ = min_ < out[i]? min_ : out[i];
    return min_;
};

template <typename type>
type vector_ocl<type>::max()
{
    // Pointer to work
    pointer input = this->copy();
    pointer output = input;

    // Working groups
    auto workGroupSize = cl_manager.get_work_group_size();
    auto numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    // std::cout << "\nwork group size =" << workGroupSize << std::endl;
    // std::cout << "num groups =" << numWorkGroups << std::endl;

    // Setup the kernel
    std::string str_kernel = kernel_minmax( string_type<type>(), true );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_max");

    while(numWorkGroups > workGroupSize)
    {
        output = minmax_loop(input, workGroupSize, numWorkGroups);
        input = output;
        numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    };

    // Finish the reduction work in CPU
    std::vector<type> out = output->std_vector();
    type max_ = out[0];
    for(int i = 0; i < out.size(); i++) max_ = max_ > out[i]? max_ : out[i];
    return max_;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::minmax_loop(typename vector_ocl<type>::pointer input, cl_int workGroupSize, cl_int numWorkGroups)
{
    auto output = vector_ocl<type>::new_pointer(numWorkGroups);

    cl_manager.get_kernel().setArg(0, *(input->get_buffer()));
    cl_manager.get_kernel().setArg(1, input->size());
    cl_manager.get_kernel().setArg(2, sizeof(type)*workGroupSize, nullptr);
    cl_manager.get_kernel().setArg(3, *(output->get_buffer()));

    cl_manager.execute(numWorkGroups*workGroupSize, workGroupSize);
    return output;
};

template <typename type>
type vector_ocl<type>::sum()
{
    // Pointer to work
    pointer input = this->copy();
    pointer output = input;

    // Working groups
    auto workGroupSize = cl_manager.get_work_group_size();
    auto numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    // std::cout << "\nwork group size =" << workGroupSize << std::endl;
    // std::cout << "num groups =" << numWorkGroups << std::endl;

    // Setup the kernel
    std::string str_kernel = kernel_sum( string_type<type>());
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_sum");

    while(numWorkGroups > workGroupSize)
    {
        output = sum_loop(input, workGroupSize, numWorkGroups);
        input = output;
        numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    };

    // Finish the reduction work in CPU
    std::vector<type> out = output->std_vector();
    type sum_ = 0;
    for(int i = 0; i < out.size(); i++) sum_ += out[i];
    return sum_;
};

template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::sum_loop(typename vector_ocl<type>::pointer input, cl_int workGroupSize, cl_int numWorkGroups)
{
    auto output = vector_ocl<type>::new_pointer(numWorkGroups);

    cl_manager.get_kernel().setArg(0, *(input->get_buffer()));
    cl_manager.get_kernel().setArg(1, sizeof(type)*workGroupSize, nullptr);
    cl_manager.get_kernel().setArg(2, *(output->get_buffer()));

    cl_manager.execute(numWorkGroups*workGroupSize, workGroupSize);
    return output;
};

// Vectorial dot product. Verify the same number of elements, then product and reduce
template <typename type>
type vector_ocl<type>::dot(const vector_ocl<type> & input)
{
    assert_size(input);
    auto output = (*this)*input;
    // output->print_data();
    return output->sum();
};

// ===========================================
// Functions
// ===========================================
template <typename type>
typename vector_ocl<type>::pointer vector_ocl<type>::normalize(type min, type max)
{
    type minv = this->min();
    type maxv = this->max();

    vector_ocl<type>::pointer output = *this - minv;
    output = *(*output*((max - min)/(maxv - minv))) + min;
    return output;
};

template <typename type> template <typename type_cast>
typename vector_ocl<type_cast>::pointer vector_ocl<type>::cast()
{
    int size = this->size();
    auto output = vector_ocl<type_cast>::new_pointer(size);
    std::string str_kernel = kernel_cast( string_type<type>(), string_type<type_cast>());
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_cast");
    cl_manager.arguments(*buffer, *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

template <typename type>
void vector_ocl<type>::pad(typename vector_ocl<type>::pointer input, typename vector_ocl<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post)
{
    if(sz.size() == 2)
    {
        std::string str_kernel = kernel_pad_2d( string_type<type>(), true );
        cl_manager.program(str_kernel, "kernel_pad_2d");
        cl_manager.arguments(*(input->get_buffer()), *(output->get_buffer()),
                                pre[0], pre[1], post[0], post[1] );
        cl_manager.execute(sz);
    }
    else if (sz.size() == 3)
    {
        std::string str_kernel = kernel_pad_3d( string_type<type>(), true );
        cl_manager.program(str_kernel, "kernel_pad_3d");
        cl_manager.arguments(*(input->get_buffer()), *(output->get_buffer()),
                                pre[0], pre[1], pre[2], post[0], post[1], post[2] );
        cl_manager.execute(sz);
    }
    else ;
};

template <typename type>
void vector_ocl<type>::unpad(typename vector_ocl<type>::pointer input, typename vector_ocl<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post)
{
    if(sz.size() == 2)
    {
        std::string str_kernel = kernel_pad_2d( string_type<type>(), false );
        cl_manager.program(str_kernel, "kernel_pad_2d");
        cl_manager.arguments(*(input->get_buffer()), *(output->get_buffer()),
                                pre[0], pre[1], post[0], post[1] );
        cl_manager.execute(sz);
    }
    else if (sz.size() == 3)
    {
        std::string str_kernel = kernel_pad_3d( string_type<type>(), false );
        cl_manager.program(str_kernel, "kernel_pad_3d");
        cl_manager.arguments(*(input->get_buffer()), *(output->get_buffer()),
                                pre[0], pre[1], pre[2], post[0], post[1], post[2] );
        cl_manager.execute(sz);
    }
    else ;
};

template <typename type>
void vector_ocl<type>::affine_2d(typename vector_ocl<type>::pointer xi,
                                 typename vector_ocl<type>::pointer yi,
                                 typename vector_ocl<type>::pointer xo,
                                 typename vector_ocl<type>::pointer yo,
                                 typename vector_ocl<type>::pointer p)
{
    std::string str_kernel = kernel_affine_2d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_affine_2d");
    cl_manager.arguments(*(xi->get_buffer()), *(yi->get_buffer()), 
                         *(xo->get_buffer()), *(yo->get_buffer()),
                         *(p->get_buffer()));
    int sz = xo->size();
    cl_manager.execute(sz); // single dim NDRange
};

template <typename type>
void vector_ocl<type>::affine_3d(typename vector_ocl<type>::pointer xi,
                                 typename vector_ocl<type>::pointer yi,
                                 typename vector_ocl<type>::pointer zi,
                                 typename vector_ocl<type>::pointer xo,
                                 typename vector_ocl<type>::pointer yo,
                                 typename vector_ocl<type>::pointer zo,
                                 typename vector_ocl<type>::pointer p)
{
    std::string str_kernel = kernel_affine_3d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_affine_3d");
    cl_manager.arguments(*(xi->get_buffer()), *(yi->get_buffer()), *(zi->get_buffer()),
                         *(xo->get_buffer()), *(yo->get_buffer()), *(zo->get_buffer()),
                         *(p->get_buffer()));
    int sz = xo->size();
    cl_manager.execute(sz); // single dim NDRange
};

template <typename type>
void vector_ocl<type>::affine_sod_2d(typename vector_ocl<type>::pointer xo,
                                     typename vector_ocl<type>::pointer yo,
                                     typename vector_ocl<type>::pointer xr,
                                     typename vector_ocl<type>::pointer yr,
                                     std::vector<double> & sod)
{
    auto p = vector_ocl<double>::new_pointer(sod.size());
    p->read_ram(sod.data(),sod.size());

    std::string str_kernel = kernel_affine_sod_2d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_affine_sod_2d");
    cl_manager.arguments(*(xo->get_buffer()), *(yo->get_buffer()), 
                         *(xr->get_buffer()), *(yr->get_buffer()),
                         *(p->get_buffer()));
    int sz = xo->size();
    cl_manager.execute(sz); // single dim NDRange
};

template <typename type>
void vector_ocl<type>::affine_sod_3d(typename vector_ocl<type>::pointer xo,
                                     typename vector_ocl<type>::pointer yo,
                                     typename vector_ocl<type>::pointer zo,
                                     typename vector_ocl<type>::pointer xr,
                                     typename vector_ocl<type>::pointer yr,
                                     typename vector_ocl<type>::pointer zr,
                                     std::vector<double> & sod)
{
    auto p = vector_ocl<double>::new_pointer(sod.size());
    p->read_ram(sod.data(),sod.size());

    std::string str_kernel = kernel_affine_sod_3d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_affine_sod_3d");
    cl_manager.arguments(*(xo->get_buffer()), *(yo->get_buffer()), *(zo->get_buffer()),
                         *(xr->get_buffer()), *(yr->get_buffer()), *(zr->get_buffer()),
                         *(p->get_buffer()));
    int sz = xo->size();
    cl_manager.execute(sz); // single dim NDRange
};

template <typename type>
void vector_ocl<type>::dfield_2d(typename vector_ocl<type>::pointer xi, typename vector_ocl<type>::pointer yi,
                                 typename vector_ocl<type>::pointer x, typename vector_ocl<type>::pointer y,
                                 typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo)
{
    std::string str_kernel = kernel_dfield_2d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_dfield_2d");
    cl_manager.arguments(*(xi->get_buffer()), *(yi->get_buffer()), 
                         *(x->get_buffer()), *(y->get_buffer()), 
                         *(xo->get_buffer()), *(yo->get_buffer()));
    int sz = xo->size();
    cl_manager.execute(sz); // single dim NDRange
};

template <typename type>
void vector_ocl<type>::dfield_3d(typename vector_ocl<type>::pointer xi, typename vector_ocl<type>::pointer yi, typename vector_ocl<type>::pointer zi,
                                 typename vector_ocl<type>::pointer x, typename vector_ocl<type>::pointer y, typename vector_ocl<type>::pointer z,
                                 typename vector_ocl<type>::pointer xo, typename vector_ocl<type>::pointer yo, typename vector_ocl<type>::pointer zo)
{
    std::string str_kernel = kernel_dfield_3d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_dfield_3d");
    cl_manager.arguments(*(xi->get_buffer()), *(yi->get_buffer()), *(zi->get_buffer()),
                         *(x->get_buffer()), *(y->get_buffer()), *(z->get_buffer()),
                         *(xo->get_buffer()), *(yo->get_buffer()), *(zo->get_buffer()));
    int sz = xo->size();
    cl_manager.execute(sz); // single dim NDRange
};

template <typename type>
std::vector<typename vector_ocl<type>::pointer> vector_ocl<type>::grid2(int w, int h, std::vector<double> & sod)
{
    int size = w*h;
    auto x = vector_ocl<type>::new_pointer(size);
    auto y = vector_ocl<type>::new_pointer(size);
    
    auto p = vector_ocl<double>::new_pointer(sod.size());
    p->read_ram(sod.data(),sod.size());

    std::string str_kernel = kernel_grid_2d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_grid_2d");
    cl_manager.arguments(*(x->get_buffer()), *(y->get_buffer()), 
                        *(p->get_buffer()), w, h);
    std::vector<int> ss({w,h});
    cl_manager.execute(ss); // multidimensional NDRange
    // cl_manager.execute(w);// single dim NDRange
    
    std::vector<typename vector_ocl<type>::pointer> xy(2);
    xy[0] = x;
    xy[1] = y;
    return xy;
};

template <typename type>
std::vector<typename vector_ocl<type>::pointer> vector_ocl<type>::grid3(int w, int h, int l, std::vector<double> & sod)
{
    int size = w*h*l;
    auto x = vector_ocl<type>::new_pointer(size);
    auto y = vector_ocl<type>::new_pointer(size);
    auto z = vector_ocl<type>::new_pointer(size);

    auto p = vector_ocl<double>::new_pointer(sod.size());
    p->read_ram(sod.data(),sod.size());

    std::string str_kernel = kernel_grid_3d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_grid_3d");
    cl_manager.arguments(*(x->get_buffer()), *(y->get_buffer()), 
                        *(z->get_buffer()), *(p->get_buffer()),
                         w, h, l);
    std::vector<int> ss({w,h,l});
    cl_manager.execute(ss); // multidimensional NDRange
    // cl_manager.execute(w);// single dim NDRange

    std::vector<typename vector_ocl<type>::pointer> xyz(3);
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    return xyz;
};

template <typename type>
void vector_ocl<type>::nearest2( typename vector_ocl<type>::pointer xo, 
                                typename vector_ocl<type>::pointer yo,
                                typename vector_ocl<type>::pointer imgr,
                                typename vector_ocl<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    std::string str_kernel = kernel_nearest_interpolation_2d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_nearest_interpolation_2d");
    cl_manager.arguments(*(xo->get_buffer()), *(yo->get_buffer()), 
                        *(imgr->get_buffer()), *(imgo->get_buffer()),
                        ref_size[0], ref_size[1]);
    cl_manager.execute(out_size); // multidimensional NDRange
};

template <typename type>
void vector_ocl<type>::nearest3( typename vector_ocl<type>::pointer xo, 
                                typename vector_ocl<type>::pointer yo,
                                typename vector_ocl<type>::pointer zo,
                                typename vector_ocl<type>::pointer imgr,
                                typename vector_ocl<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    std::string str_kernel = kernel_nearest_interpolation_3d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_nearest_interpolation_3d");
    cl_manager.arguments(*(xo->get_buffer()), *(yo->get_buffer()), *(zo->get_buffer()),
                        *(imgr->get_buffer()), *(imgo->get_buffer()),
                        ref_size[0], ref_size[1], ref_size[2]);
    cl_manager.execute(out_size); // multidimensional NDRange
};

template <typename type>
void vector_ocl<type>::linear2( typename vector_ocl<type>::pointer xo, 
                                typename vector_ocl<type>::pointer yo,
                                typename vector_ocl<type>::pointer imgr,
                                typename vector_ocl<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    std::string str_kernel = kernel_linear_interpolation_2d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_linear_interpolation_2d");
    cl_manager.arguments(*(xo->get_buffer()), *(yo->get_buffer()), 
                        *(imgr->get_buffer()), *(imgo->get_buffer()),
                        ref_size[0], ref_size[1]);
    cl_manager.execute(out_size); // multidimensional NDRange
};

template <typename type>
void vector_ocl<type>::linear3( typename vector_ocl<type>::pointer xo, 
                                typename vector_ocl<type>::pointer yo,
                                typename vector_ocl<type>::pointer zo,
                                typename vector_ocl<type>::pointer imgr,
                                typename vector_ocl<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    std::string str_kernel = kernel_linear_interpolation_3d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_linear_interpolation_3d");
    cl_manager.arguments(*(xo->get_buffer()), *(yo->get_buffer()), *(zo->get_buffer()),
                        *(imgr->get_buffer()), *(imgo->get_buffer()),
                        ref_size[0], ref_size[1], ref_size[2]);
    cl_manager.execute(out_size); // multidimensional NDRange
};

template <typename type>
void vector_ocl<type>::cubic2(  typename vector_ocl<type>::pointer xo, 
                                typename vector_ocl<type>::pointer yo,
                                typename vector_ocl<type>::pointer imgr,
                                typename vector_ocl<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    std::string str_kernel = kernel_cubic_interpolation_2d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_cubic_interpolation_2d");
    cl_manager.arguments(*(xo->get_buffer()), *(yo->get_buffer()), 
                        *(imgr->get_buffer()), *(imgo->get_buffer()),
                        ref_size[0], ref_size[1]);
    cl_manager.execute(out_size); // multidimensional NDRange
};

template <typename type>
void vector_ocl<type>::cubic3( typename vector_ocl<type>::pointer xo, 
                                typename vector_ocl<type>::pointer yo,
                                typename vector_ocl<type>::pointer zo,
                                typename vector_ocl<type>::pointer imgr,
                                typename vector_ocl<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    std::string str_kernel = kernel_cubic_interpolation_3d( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_cubic_interpolation_3d");
    cl_manager.arguments(*(xo->get_buffer()), *(yo->get_buffer()), *(zo->get_buffer()),
                        *(imgr->get_buffer()), *(imgo->get_buffer()),
                        ref_size[0], ref_size[1], ref_size[2]);
    cl_manager.execute(out_size); // multidimensional NDRange
};

template <typename type>
void vector_ocl<type>::fft(std::vector<pointer> & input, std::vector<pointer> & output, std::vector<int> size, bool forward)
{
    // timer t("ms");
    // t.start();
    // Dimension
    cl_int err;
    int dim = size.size();
    size_t lengths[dim];

    // Input and Output buffer
    cl_mem buffers_in[2];
    cl_mem buffers_out[2];
    // input[0]->print_data("vin0 ocl");
    // input[1]->print_data("vin1 ocl");

    // Initiliaze
    for(int i = 0; i < dim; i++) lengths[i] = size[i];
    for(int i = 0; i < 2; i++) // real and img
    {
        buffers_in[i] = (*(input[i]->get_buffer()))();
        buffers_out[i] = (*(output[i]->get_buffer()))();
        // std::cout << "in buffer " << i << ": " << (*(input[i]->get_buffer()))() << std::endl;
        // std::cout << "out buffer " << i << ": " << (*(output[i]->get_buffer()))() << std::endl;

    };
    // t.lap("[Init]");

    // Temporary buffer
    cl_mem tmp_buffer = 0;
    size_t tmp_buffer_size = 0;
    int status = 0;
    
    // Setup clFFT
    clfftSetupData fft_setup;
    err = clfftInitSetupData(&fft_setup);
    err = clfftSetup(&fft_setup);
    // t.lap("[Setup]");

    // Create a default plan for a complex FFT
    clfftPlanHandle plan_handle;
    clfftDim dd;
    if (dim == 2) dd = CLFFT_2D;
    if (dim == 3) dd = CLFFT_3D;

    // t.lap("[Handle]");
    cl_context ctx = (cl_manager.get_context())();
    err = clfftCreateDefaultPlan(&plan_handle, ctx, dd, lengths);
    // t.lap("[Plan]");

    // Set plan parameters
    if (string_type<type>() == "double") err = clfftSetPlanPrecision(plan_handle, CLFFT_DOUBLE);
    else err = clfftSetPlanPrecision(plan_handle, CLFFT_SINGLE);    
    err = clfftSetLayout(plan_handle, CLFFT_COMPLEX_PLANAR, CLFFT_COMPLEX_PLANAR);
    err = clfftSetResultLocation(plan_handle, CLFFT_OUTOFPLACE);
    // t.lap("[PlanProp]");

    // Bake the plan
    cl_command_queue queue = (cl_manager.get_queue())();
    err = clfftBakePlan(plan_handle, 1, &queue, NULL, NULL);
    // t.lap("[Bake]");

    // Create temporary buffer
    status = clfftGetTmpBufSize(plan_handle, &tmp_buffer_size);

    if ((status == 0) && (tmp_buffer_size > 0))
    {
        // tmp_buffer = cl::Buffer(cl_manager.get_context(), CL_MEM_READ_WRITE, tmp_buffer_size, nullptr, &err);
        tmp_buffer = clCreateBuffer(ctx, CL_MEM_READ_WRITE, tmp_buffer_size, 0, &err);
        if (err != CL_SUCCESS)
            printf("Error with tmp_buffer clCreateBuffer\n");
    };
    // t.lap("[Buffer]");
    
    // Execute Forward or Backward FFT
    if(forward)
    {
        // Execute the plan
        err = clfftEnqueueTransform(plan_handle, CLFFT_FORWARD, 1, &queue, 0, NULL, NULL,
        buffers_in, buffers_out, tmp_buffer);
        if (err != CL_SUCCESS)
            printf("Error running fft plan\n");
    }
    else
    {
        // Execute the plan
        err = clfftEnqueueTransform(plan_handle, CLFFT_BACKWARD, 1, &queue, 0, NULL, NULL,
        buffers_in, buffers_out, tmp_buffer);
        if (err != CL_SUCCESS)
            printf("Error running fft plan\n");
    };
    // t.lap("[FFT]");
    // t.finish();

    // Release the plan
    err = clfftDestroyPlan(&plan_handle);

    // Release clFFT library
    clfftTeardown();
};

template <typename type>
void vector_ocl<type>::gradientx( typename vector_ocl<type>::pointer imgr,
                                  typename vector_ocl<type>::pointer imgo,
                                  std::vector<int> ref_size)
{
    std::string str_kernel = kernel_gradientx( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_gradientx");
    cl_manager.arguments(*(imgr->get_buffer()), *(imgo->get_buffer()));
    cl_manager.execute(ref_size); // multidimensional NDRange
};

template <typename type>
void vector_ocl<type>::gradienty( typename vector_ocl<type>::pointer imgr,
                                  typename vector_ocl<type>::pointer imgo,
                                  std::vector<int> ref_size)
{
    std::string str_kernel = kernel_gradienty( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_gradienty");
    cl_manager.arguments(*(imgr->get_buffer()), *(imgo->get_buffer()));
    cl_manager.execute(ref_size); // multidimensional NDRange
};

template <typename type>
void vector_ocl<type>::gradientz( typename vector_ocl<type>::pointer imgr,
                                  typename vector_ocl<type>::pointer imgo,
                                  std::vector<int> ref_size)
{
    std::string str_kernel = kernel_gradientz( string_type<type>() );
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_gradientz");
    cl_manager.arguments(*(imgr->get_buffer()), *(imgo->get_buffer()));
    cl_manager.execute(ref_size); // multidimensional NDRange
};

template <typename type>
void vector_ocl<type>::convolution( typename vector_ocl<type>::pointer imgr,
                                    typename vector_ocl<type>::pointer kernel,
                                    typename vector_ocl<type>::pointer imgo,
                                    std::vector<int> ref_size, int kwidth)
{
    if (ref_size.size() == 2)
    {
        std::string str_kernel = kernel_convolution_2d( string_type<type>() );
        // std::cout << str_kernel << std::endl;
        cl_manager.program(str_kernel, "kernel_convolution_2d");
        cl_manager.arguments(*(imgr->get_buffer()), *(kernel->get_buffer()), 
                            *(imgo->get_buffer()), kwidth);
        cl_manager.execute(ref_size); // multidimensional NDRange
    }
    else if (ref_size.size() == 3)
    {
        std::string str_kernel = kernel_convolution_3d( string_type<type>() );
        // std::cout << str_kernel << std::endl;
        cl_manager.program(str_kernel, "kernel_convolution_3d");
        cl_manager.arguments(*(imgr->get_buffer()), *(kernel->get_buffer()), 
                            *(imgo->get_buffer()), kwidth);
        cl_manager.execute(ref_size); // multidimensional NDRange
    }
    else ;
};

}; //end namespace

#endif