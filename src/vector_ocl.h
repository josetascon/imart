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
#include <cassert>      // assert

// gpu libs
#include <CL/cl.hpp>

// local libs
#include "object.h"
#include "kernels.h"
#include "opencl_object.h"

namespace imart
{

opencl_object cl_manager; // manager of vector ocl;

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
    int size_;
    cl_int err;
    std::shared_ptr<cl::Buffer> buffer;

    // ===========================================
    // Constructor Functions
    // ===========================================
    void init(int s);
    void allocate(int s);
    
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
    void read_ram(type * p, int size);
    void write_ram(type * p, int size);

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

    // // ===========================================
    // // Reduction functions
    // // ===========================================
    type min();
    type max();
    type sum();
    // // type prod();     // may produce overflow error
    // type dot(const vector_ocl<type> & input);

    // // ===========================================
    // // Functions
    // // ===========================================
    // pointer normalize(type min = 0.0, type max = 1.0);

    template<typename type_cast>
    typename vector_ocl<type_cast>::pointer cast();

    static std::vector<typename vector_ocl<type>::pointer> grid_2d(int w, int h, std::vector<double> & sod);
    static std::vector<typename vector_ocl<type>::pointer> grid_3d(int w, int h, int l, std::vector<double> & sod);
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
    size_ = s;
    if(size_>0) allocate(size_);
    else buffer = nullptr;
};

template <typename type>
void vector_ocl<type>::allocate(int s)
{
    err = 0;
    buffer = std::make_shared<cl::Buffer>(cl_manager.get_context(), CL_MEM_READ_WRITE, sizeof(type)*s, nullptr, &err);
    assert(err == 0);
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
    cl_manager.program(str_kernel, "kernel_copy");
    cl_manager.arguments(*(input.get_buffer()), *buffer);
    cl_manager.execute(size_);
};

template <typename type>
void vector_ocl<type>::copy_(const vector_ocl<type> & input)
{
    // std::cout << "copy_";
    size_ = input.size();
    buffer.reset();
    buffer = input.get_buffer();
};

template <typename type>
void vector_ocl<type>::mimic_(const vector_ocl<type> & input)
{
    // std::cout << "mimic_";
    size_ = input.size();
    buffer.reset();
    allocate(size_);
};

// ===========================================
// Get Functions
// ===========================================
template <typename type>
int vector_ocl<type>::size() const
{
    return size_;
};

template <typename type>
std::shared_ptr<cl::Buffer> vector_ocl<type>::get_buffer() const
{
    return buffer;
};

template <typename type>
std::vector<type> vector_ocl<type>::std_vector()
{
    std::vector<type> tmp(size_);
    if (size_ > 0)
        cl_manager.get_queue().enqueueReadBuffer(*buffer, CL_TRUE, 0, sizeof(type)*tmp.size(), tmp.data());
    return tmp;
};

// ===========================================
// Memory Functions
// ===========================================
template <typename type>
void vector_ocl<type>::read_ram(type * p, int s)
{
    if(s > 0)
        cl_manager.get_queue().enqueueWriteBuffer(*buffer, CL_TRUE, 0, sizeof(type)*(s), p);
};

template <typename type>
void vector_ocl<type>::write_ram(type * p, int s)
{
    if(s > 0)
        cl_manager.get_queue().enqueueReadBuffer(*buffer, CL_TRUE, 0, sizeof(type)*s, p);
};


// ===========================================
// Print Functions
// ===========================================
template <typename type>
std::string vector_ocl<type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Vector GPU Information";
    if (msg != "") { title = msg; };

    // Summary of the object information
    ss << object::info(title);
    ss << "Size: \t\t\t" << size_ << std::endl;
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
    cl_manager.execute(size_);
};

template <typename type>
void vector_ocl<type>::random(float min, float max)
{
    std::string str_kernel = kernel_random( string_type<type>(), min, max);
    // // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_random");
    cl_manager.arguments(*buffer);
    cl_manager.execute(size_);
};

// ===========================================
// Overloading operators
// ===========================================
template <typename type>
type vector_ocl<type>::operator [](int e)
{
    std::vector<type> tmp(1);
    if(e >= 0 && e < size_ && size_ > 0)
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
    cl_manager.execute(size_);
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
    cl_manager.execute(size_);
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
    cl_manager.execute(size_);
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
    cl_manager.execute(size_);
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
    cl_manager.execute(size_);
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
    cl_manager.execute(size_);
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
    cl_manager.execute(size_);
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
    cl_manager.execute(size_);
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
    cl_manager.execute(size_);
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
    cl_manager.execute(size_);
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
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "-", false, true);  // function = true, reverse = true
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
    
    std::string str_kernel = kernel_scalar( string_type<type>(), "/", false, true);  // function = true, reverse = true
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(*(input.get_buffer()), *(output->get_buffer()), scalar);
    cl_manager.execute(size);
    return output;
};

// ===========================================
// Reduction Functions
// ===========================================
template <typename type>
type vector_ocl<type>::min()
{
    std::string str_kernel = kernel_minmax( string_type<type>(), false);
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_min");
    
    auto workGroupSize = cl_manager.get_work_group_size();
    auto numWorkGroups = 1 + ((size_ - 1)/workGroupSize);
    // std::cout << "\nwork group size =" << workGroupSize << std::endl;
    // std::cout << "num groups =" << numWorkGroups << std::endl;
    
    auto output = vector_ocl<type>::new_pointer(numWorkGroups);

    cl_manager.get_kernel().setArg(0, *buffer);
    cl_manager.get_kernel().setArg(1, size_);
    cl_manager.get_kernel().setArg(2, sizeof(type)*workGroupSize, nullptr);
    cl_manager.get_kernel().setArg(3, *(output->get_buffer()));

    cl_manager.get_queue().enqueueNDRangeKernel(cl_manager.get_kernel(), cl::NullRange, cl::NDRange(numWorkGroups*workGroupSize), cl::NDRange(workGroupSize));

    std::vector<type> out = output->std_vector();

    // type * p = out.data();
    // std::cout << "[ ";
    // for(int k=0; k<out.size(); k++) { std::cout << p[k] << " "; };
    // std::cout << "]" << std::endl;

    // Finish the reduction work in CPU
    type min_ = out[0];
    for(int i = 0; i < out.size(); i++) min_ = min_ < out[i]? min_ : out[i];
    return min_;
};

template <typename type>
type vector_ocl<type>::max()
{
    std::string str_kernel = kernel_minmax( string_type<type>(), true);
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_max");
    
    auto workGroupSize = cl_manager.get_work_group_size();
    auto numWorkGroups = 1 + ((size_ - 1)/workGroupSize);
    // std::cout << "\nwork group size =" << workGroupSize << std::endl;
    // std::cout << "num groups =" << numWorkGroups << std::endl;
    
    auto output = vector_ocl<type>::new_pointer(numWorkGroups);

    cl_manager.get_kernel().setArg(0, *buffer);
    cl_manager.get_kernel().setArg(1, size_);
    cl_manager.get_kernel().setArg(2, sizeof(type)*workGroupSize, nullptr);
    cl_manager.get_kernel().setArg(3, *(output->get_buffer()));

    cl_manager.get_queue().enqueueNDRangeKernel(cl_manager.get_kernel(), cl::NullRange, cl::NDRange(numWorkGroups*workGroupSize), cl::NDRange(workGroupSize));

    std::vector<type> out = output->std_vector();

    // type * p = out.data();
    // std::cout << "[ ";
    // for(int k=0; k<out.size(); k++) { std::cout << p[k] << " "; };
    // std::cout << "]" << std::endl;

    // Finish the reduction work in CPU
    type max_ = out[0];
    for(int i = 0; i < out.size(); i++) max_ = max_ > out[i]? max_ : out[i];
    return max_;
};

template <typename type>
type vector_ocl<type>::sum()
{
    std::string str_kernel = kernel_sum( string_type<type>());
    // std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_sum");
    
    auto workGroupSize = cl_manager.get_work_group_size();
    auto numWorkGroups = 1 + ((size_ - 1)/workGroupSize);
    // std::cout << "\nwork group size =" << workGroupSize << std::endl;
    // std::cout << "num groups =" << numWorkGroups << std::endl;

    auto output = vector_ocl<type>::new_pointer(numWorkGroups);

    cl_manager.get_kernel().setArg(0, *buffer);
    cl_manager.get_kernel().setArg(1, sizeof(type)*workGroupSize, nullptr);
    cl_manager.get_kernel().setArg(2, *(output->get_buffer()));

    cl_manager.get_queue().enqueueNDRangeKernel(cl_manager.get_kernel(), cl::NullRange, cl::NDRange(numWorkGroups*workGroupSize), cl::NDRange(workGroupSize));

    std::vector<type> out = output->std_vector();

    // type * p = out.data();
    // std::cout << "[ ";
    // for(int k=0; k<out.size(); k++) { std::cout << p[k] << " "; };
    // std::cout << "]" << std::endl;

    // Finish the reduction work in CPU
    type sum_ = 0;
    for(int i = 0; i < out.size(); i++) sum_ += out[i];
    return sum_;
};

// Vectorial dot product. Verify the same number of elements, then product and reduce
// template <typename type>
// type vector_ocl<type>::dot(const vector_ocl<type> & input)
// {
//     assert_size(input);

//     type x;
//     viennacl::scalar<type> xm(0);
//     viennacl::linalg::inner_prod_impl(*this, input, xm);
//     x = xm;
//     return x;
// };

// ===========================================
// Functions
// ===========================================
// template <typename type>
// typename vector_ocl<type>::pointer vector_ocl<type>::normalize(type min, type max)
// {
//     type minv = this->min();
//     type maxv = this->max();

//     vector_ocl<type>::pointer output = *this - minv;
//     output = *(*output*((max - min)/(maxv - minv))) + min;
//     return output->clone(); // memory error if viennacl is out of context, clone for keeping alive!!!
// };

template <typename type> template <typename type_cast>
typename vector_ocl<type_cast>::pointer vector_ocl<type>::cast()
{
    int size = this->size();
    auto output = vector_ocl<type_cast>::new_pointer(size);
    // viennacl::linalg::convert(*output,*this);
    std::string str_kernel = kernel_cast( string_type<type>(), string_type<type_cast>());
    cl_manager.program(str_kernel, "kernel_cast");
    cl_manager.arguments(*buffer, *(output->get_buffer()));
    cl_manager.execute(size);
    return output;
};

// template <typename type>
// std::vector<typename vector_ocl<type>::pointer> grid_2d(int w, int h)
template <typename type>
std::vector<typename vector_ocl<type>::pointer> vector_ocl<type>::grid_2d(int w, int h, std::vector<double> & sod)
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
    cl_manager.execute(size);
    
    std::vector<typename vector_ocl<type>::pointer> xy(2);
    xy[0] = x;
    xy[1] = y;
    return xy;
};

template <typename type>
std::vector<typename vector_ocl<type>::pointer> vector_ocl<type>::grid_3d(int w, int h, int l, std::vector<double> & sod)
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
    cl_manager.execute(size);
    
    std::vector<typename vector_ocl<type>::pointer> xyz(3);
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    return xyz;
};

}; //end namespace

#endif