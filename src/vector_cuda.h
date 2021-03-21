/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __VECTOR_CUDA_H__
#define __VECTOR_CUDA_H__

// std libs
#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <random>       // std::random
#include <cassert>      // assert

// local libs
#include "object.h"
#include "cuda_object.h"
#include "cuda_buffer.h"
#include "utils/timer.h"
#include "vector_cpu.h"

namespace imart
{

// #ifndef IMART_WITH_CUDA
// void cuda_check_gpu()
// {
//     std::cout << "CUDA is not available with current imart built" << std::endl;
// }
// #endif

#ifdef IMART_WITH_CUDA
cuda_object cuda_manager; // manager of vector cuda;

// imart function to check cuda device
void imart_cuda_device_name()
{
    cuda_check_gpu(); // TODO: find a way to leave it when cuda is disabled
};
#endif

// Class object
template <typename type>
class vector_cuda: public inherit<vector_cuda<type>, object>
{
public:
    //Type definitions
    using self    = vector_cuda;
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
    int err;
    std::shared_ptr<cuda_buffer<type>> buffer;
    vector_cpu<type> data_cpu;

    // ===========================================
    // Constructor Functions
    // ===========================================
    void init(int s);
    void allocate(int s);

    // ===========================================
    // Reduction functions
    // ===========================================
    pointer sum_loop(pointer input, int workGroupSize, int numWorkGroups);
    pointer minmax_loop(pointer input, int workGroupSize, int numWorkGroups, bool is_max);
    
public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    vector_cuda();                               // constructor empty
    vector_cuda(int s);                          // constructor size
    vector_cuda(int s, type value);              // constructor assign
    vector_cuda(std::initializer_list<type> list);
    vector_cuda(const vector_cuda & input);       // constructor clone
    ~vector_cuda();                              // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const vector_cuda & input);  // copy everything
    virtual void copy_ (const vector_cuda & input);  // share data
    virtual void mimic_(const vector_cuda & input);  // copy meta data

    // ===========================================
    // Get Functions
    // ===========================================
    int size() const;
    type * data() const;
    std::shared_ptr<cuda_buffer<type>> get_buffer() const;
    std::vector<type> std_vector();

    // ===========================================
    // Memory Functions
    // ===========================================
    void read_ram(type * p, int size, int offset = 0);
    void write_ram(type * p, int size, int offset = 0);
    void equal(pointer input);
    void to_cpu();

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    void assert_size(const vector_cuda<type> & input);

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
    vector_cuda<type> & operator = (const vector_cuda & input);
    pointer operator = (const vector_cuda<type>::pointer input);
    pointer operator + (const vector_cuda<type> & input);
    pointer operator - (const vector_cuda<type> & input);
    pointer operator * (const vector_cuda<type> & input);
    pointer operator / (const vector_cuda<type> & input);
    pointer operator ^ (const vector_cuda<type> & input);

    // // Scalar to vector
    pointer operator + (type scalar);
    pointer operator - (type scalar);
    pointer operator * (type scalar);
    pointer operator / (type scalar);
    pointer operator ^ (type scalar);

    // Friend classes to support scalar to vector left hand side
    template<typename type_>
    friend typename vector_cuda<type_>::pointer operator + (type_ scalar, vector_cuda<type_> & input);
    template<typename type_>
    friend typename vector_cuda<type_>::pointer operator - (type_ scalar, const vector_cuda<type_> & input);
    template<typename type_>
    friend typename vector_cuda<type_>::pointer operator * (type_ scalar, vector_cuda<type_> & input);
    template<typename type_>
    friend typename vector_cuda<type_>::pointer operator / (type_ scalar, const vector_cuda<type_> & input);
    template<typename type_>
    friend typename vector_cuda<type_>::pointer operator ^ (type_ scalar, const vector_cuda<type_> & input);

    // Conditions 
    pointer operator == (const vector_cuda<type> & input);
    pointer operator > (const vector_cuda<type> & input);
    pointer operator < (const vector_cuda<type> & input);
    pointer operator >= (const vector_cuda<type> & input);
    pointer operator <= (const vector_cuda<type> & input);

    pointer operator == (type scalar);
    pointer operator > (type scalar);
    pointer operator < (type scalar); // std::shared_ptr<vector_cuda<type>>
    pointer operator >= (type scalar);
    pointer operator <= (type scalar);

    void replace(const vector_cuda<type> & idxs, const vector_cuda<type> & input);
    void replace(const vector_cuda<type> & idxs, type value);

    // ===========================================
    // Reduction functions
    // ===========================================
    type min();
    type max();
    type sum();
    type dot(const vector_cuda<type> & input);
    // type prod();     // may produce overflow error

    // // ===========================================
    // // Functions
    // // ===========================================
    pointer normalize(type min = 0.0, type max = 1.0);

    template<typename type_cast>
    typename vector_cuda<type_cast>::pointer cast();

    static void pad(typename vector_cuda<type>::pointer input, typename vector_cuda<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post);
    static void unpad(typename vector_cuda<type>::pointer input, typename vector_cuda<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post);

    static void affine_2d(typename vector_cuda<type>::pointer xi, typename vector_cuda<type>::pointer yi,
                          typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo,
                          typename vector_cuda<type>::pointer p);

    static void affine_3d(typename vector_cuda<type>::pointer xi, typename vector_cuda<type>::pointer yi, typename vector_cuda<type>::pointer zi,
                          typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo, typename vector_cuda<type>::pointer zo,
                          typename vector_cuda<type>::pointer p);

    static void affine_sod_2d(typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo,
                              typename vector_cuda<type>::pointer xr, typename vector_cuda<type>::pointer yr, 
                              std::vector<double> & sod);
    static void affine_sod_3d(typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo, typename vector_cuda<type>::pointer zo,
                              typename vector_cuda<type>::pointer xr, typename vector_cuda<type>::pointer yr, typename vector_cuda<type>::pointer zr,
                              std::vector<double> & sod);

    static void dfield_2d(typename vector_cuda<type>::pointer xi, typename vector_cuda<type>::pointer yi,
                          typename vector_cuda<type>::pointer x, typename vector_cuda<type>::pointer y,
                          typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo);

    static void dfield_3d(typename vector_cuda<type>::pointer xi, typename vector_cuda<type>::pointer yi, typename vector_cuda<type>::pointer zi,
                          typename vector_cuda<type>::pointer x, typename vector_cuda<type>::pointer y, typename vector_cuda<type>::pointer z,
                          typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo, typename vector_cuda<type>::pointer zo);
    
    static std::vector<typename vector_cuda<type>::pointer> grid2(int w, int h, std::vector<double> & sod);
    static std::vector<typename vector_cuda<type>::pointer> grid3(int w, int h, int l, std::vector<double> & sod);

    static void nearest2( typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo,
                        typename vector_cuda<type>::pointer imgr, typename vector_cuda<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);
    static void nearest3( typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo, typename vector_cuda<type>::pointer zo,
                        typename vector_cuda<type>::pointer imgr, typename vector_cuda<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);

    static void linear2( typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo,
                        typename vector_cuda<type>::pointer imgr, typename vector_cuda<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);
    static void linear3( typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo, typename vector_cuda<type>::pointer zo,
                        typename vector_cuda<type>::pointer imgr, typename vector_cuda<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);

    static void cubic2( typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo,
                        typename vector_cuda<type>::pointer imgr, typename vector_cuda<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);
    static void cubic3( typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo, typename vector_cuda<type>::pointer zo,
                        typename vector_cuda<type>::pointer imgr, typename vector_cuda<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);

    static void fft(std::vector<pointer> & input, std::vector<pointer> & output, std::vector<int> size, bool forward);

    static void gradientx( typename vector_cuda<type>::pointer imgr,
                           typename vector_cuda<type>::pointer imgo,
                           std::vector<int> ref_size);

    static void gradienty( typename vector_cuda<type>::pointer imgr,
                           typename vector_cuda<type>::pointer imgo,
                           std::vector<int> ref_size);

    static void gradientz( typename vector_cuda<type>::pointer imgr,
                           typename vector_cuda<type>::pointer imgo,
                           std::vector<int> ref_size);

    static void convolution( typename vector_cuda<type>::pointer imgr,
                             typename vector_cuda<type>::pointer kernel,
                             typename vector_cuda<type>::pointer imgo,
                             std::vector<int> ref_size, std::vector<int> kernel_size);
};


// ===========================================
//          Functions of Class vector_cuda
// ===========================================

// ===========================================
// Constructors
// ===========================================
template <typename type>
vector_cuda<type>::vector_cuda()
{
    class_name = "vector_cuda";
    init(0);
};

template <typename type>
vector_cuda<type>::vector_cuda(int s)
{
    class_name = "vector_cuda";
    init(s);
};

template <typename type>
vector_cuda<type>::vector_cuda(int s, type value)
{
    class_name = "vector_cuda";
    init(s);
    assign(value);
};

template <typename type>
vector_cuda<type>::vector_cuda(std::initializer_list<type> list)
{
    // printf("constructor init list\n");
    class_name = "vector_cuda";
    int s = list.size();
    init(s);
    if (s > 0)
        buffer->push_memory(list.begin(),s);
        // cuda_manager.get_queue().enqueueWriteBuffer(*buffer, CL_TRUE, 0, sizeof(type)*(s), list.begin());
};

template <typename type>
vector_cuda<type>::vector_cuda(const vector_cuda<type> & input)
{
    class_name = "vector_cuda";
    clone_(input);          // call the virtual
};

//! Destructor
template <typename type>
vector_cuda<type>::~vector_cuda()
{
    ;
};

template <typename type>
void vector_cuda<type>::init(int s)
{
    _size_ = s;
    if(_size_>0) allocate(s);
    else buffer = nullptr;
};

template <typename type>
void vector_cuda<type>::allocate(int s)
{
    // err = 0;
    buffer = std::make_shared<cuda_buffer<type>>(s);
    data_cpu = vector_cpu<type>(s);
    // assert(err == 0);
};

// ===========================================
// Create Functions
// ===========================================
template <typename type>
void vector_cuda<type>::clone_(const vector_cuda<type> & input)
{
    // std::cout << "clone_";
    mimic_(input);
    // clone buffer;
    cuda_manager.setup(_size_);
    cuda_manager.execute(cuda_kernel_copy<type>, input.get_buffer()->get(), buffer->get(), _size_);
};

template <typename type>
void vector_cuda<type>::copy_(const vector_cuda<type> & input)
{
    // std::cout << "copy_";
    _size_ = input.size();
    buffer.reset();
    buffer = input.get_buffer();
};

template <typename type>
void vector_cuda<type>::mimic_(const vector_cuda<type> & input)
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
int vector_cuda<type>::size() const
{
    return _size_;
};

template <typename type>
type * vector_cuda<type>::data() const
{
    return data_cpu.data();
};

template <typename type>
std::shared_ptr<cuda_buffer<type>> vector_cuda<type>::get_buffer() const
{
    return buffer;
};

template <typename type>
std::vector<type> vector_cuda<type>::std_vector()
{
    std::vector<type> tmp(_size_);
    if (_size_ > 0)
        buffer->pull_memory(tmp.data(), _size_);
        // std::cout << "std_vector" << std::endl;
    return tmp;
};

// ===========================================
// Memory Functions
// ===========================================
template <typename type>
void vector_cuda<type>::read_ram(type * p, int s, int offset)
{
    if(s > 0)
        buffer->push_memory(p, s, offset);
};

template <typename type>
void vector_cuda<type>::write_ram(type * p, int s, int offset)
{
    if(s > 0)
        buffer->pull_memory(p, s, offset);
};

template <typename type>
void vector_cuda<type>::equal(pointer input)
{
    // Used to hold a cpu pointer. In gpu there is no need to hold the pointer.
    // Therefore, we use this to write to cpu
    
    assert_size(*input);

    // GPU
    buffer.reset();
    buffer = input->get_buffer();

    // Tranfer GPU to CPU (keeping same CPU pointer)
    to_cpu();
};

template <typename type>
void vector_cuda<type>::to_cpu()
{
    // Tranfer GPU to CPU (keeping same CPU pointer)
    write_ram( data_cpu.data(), _size_ );
};

// ===========================================
// Print Functions
// ===========================================
template <typename type>
std::string vector_cuda<type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Vector CUDA Information";
    if (msg != "") { title = msg; };

    // Summary of the object information
    ss << object::info(title);
    ss << "Size: \t\t\t" << _size_ << std::endl;
    ss << "Memory: \t\t" << buffer << std::endl;
    if (buffer != nullptr) ss << "Address cuda: \t\t" << buffer->get() << std::endl;
    else ss << "Address cuda: \t\t" << 0 << std::endl;
    // ss << "Capacity: \t\t" << this->capacity() << std::endl;
    return ss.str();
};

template <typename type>
std::string vector_cuda<type>::info_data(std::string msg)
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
void vector_cuda<type>::assert_size(const vector_cuda<type> & input)
{
    imart_assert(this->size() == input.size(), "Mismatch of vector size");
};

// ===========================================
// Initialization Functions
// ===========================================
template <typename type>
void vector_cuda<type>::zeros()
{
    assign(type(0));
};

template <typename type>
void vector_cuda<type>::ones()
{
    assign(type(1));
};

// Set all pixel to a fixed value
template <typename type>
void vector_cuda<type>::assign(type value)
{
    if( _size_ > 0 )
    {
        type v = value;
        cuda_manager.setup(_size_);
        cuda_manager.execute(cuda_kernel_assign<type>, buffer->get(), v, _size_);
    };

    // manual call
    // auto block = std::vector<int>{32,0,0};
    // auto grid = std::vector<int>{1,0,0};
    // cuda_kernel_assign<type>(grid,block,buffer->get(), value, _size_);
};

template <typename type>
void vector_cuda<type>::random(float min, float max)
{
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
type vector_cuda<type>::operator [](int e)
{
    std::vector<type> tmp(1);
    if(e >= 0 && e < _size_ && _size_ > 0)
        // std::cout << "cl read" << std::endl;
        buffer->pull_memory(tmp.data(), 1, e); //offset e
    return tmp[0];
};

template <typename type>
vector_cuda<type> & vector_cuda<type>::operator = (const vector_cuda & input)
{
    copy_(input);
    return *this;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator = (const vector_cuda<type>::pointer input)
{
    return input;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator + (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);

    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_add<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator - (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);

    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_sub<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator * (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);

    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_mul<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator / (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);

    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_div<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator ^ (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_pow<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

// Scalar right hand side
template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator + (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);

    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_add_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator - (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_sub_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator * (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_mul_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator / (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_div_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator ^ (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_pow_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

// Scalar left hand side
template <typename type>
typename vector_cuda<type>::pointer operator + (type scalar, vector_cuda<type> & input)
{
    return input + scalar;
};

template <typename type>
typename vector_cuda<type>::pointer operator - (type scalar, const vector_cuda<type> & input)
{
    int size = input.size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_sub_scalar_inv<type>, input.get_buffer()->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer operator * (type scalar, vector_cuda<type> & input)
{
    return input * scalar;
};

template <typename type>
typename vector_cuda<type>::pointer operator / (type scalar, const vector_cuda<type> & input)
{
    int size = input.size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_div_scalar_inv<type>, input.get_buffer()->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer operator ^ (type scalar, const vector_cuda<type> & input)
{
    int size = input.size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_pow_scalar_inv<type>, input.get_buffer()->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator == (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_equal<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator > (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_greater<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator < (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_less<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator >= (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_greater_equal<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator <= (const vector_cuda<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_less_equal<type>, buffer->get(), input.get_buffer()->get(), output->get_buffer()->get(), size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator == (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_equal_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator > (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_greater_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator < (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_less_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator >= (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_greater_equal_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::operator <= (type scalar)
{
    int size = this->size();
    auto output = vector_cuda<type>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_less_equal_scalar<type>, buffer->get(), output->get_buffer()->get(), scalar, size);
    return output;
};

template <typename type>
void vector_cuda<type>::replace(const vector_cuda<type> & idxs, const vector_cuda<type> & input)
{
    assert_size(idxs);
    assert_size(input);
    int size = this->size();

    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_replace<type>, idxs.get_buffer()->get(), input.get_buffer()->get(), buffer->get(), size );
    return;
};

template <typename type>
void vector_cuda<type>::replace(const vector_cuda<type> & idxs, type value)
{
    assert_size(idxs);
    int size = this->size();

    cuda_manager.setup(size);
    cuda_manager.execute(cuda_kernel_replace_scalar<type>, idxs.get_buffer()->get(), buffer->get(), value, size );
    return;
};

// ===========================================
// Reduction Functions
// ===========================================
template <typename type>
type vector_cuda<type>::min()
{
    // Pointer to work
    pointer input = this->copy();
    pointer output = input;

    // Working groups
    auto workGroupSize = cuda_manager.get_threads();
    auto numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    // std::cout << "\nwork group size =" << workGroupSize << std::endl;
    // std::cout << "num groups =" << numWorkGroups << std::endl;

    // Main loop
    while(numWorkGroups > workGroupSize)
    {
        output = minmax_loop(input, workGroupSize, numWorkGroups, false);
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
type vector_cuda<type>::max()
{
    // Pointer to work
    pointer input = this->copy();
    pointer output = input;

    // Working groups
    auto workGroupSize = cuda_manager.get_threads();
    auto numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    // std::cout << "\nwork group size =" << workGroupSize << std::endl;
    // std::cout << "num groups =" << numWorkGroups << std::endl;

    // Main loop
    while(numWorkGroups > workGroupSize)
    {
        output = minmax_loop(input, workGroupSize, numWorkGroups, true);
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
typename vector_cuda<type>::pointer vector_cuda<type>::minmax_loop(typename vector_cuda<type>::pointer input, int workGroupSize, int numWorkGroups, bool is_max)
{
    auto output = vector_cuda<type>::new_pointer(numWorkGroups);

    if (is_max)
    {
        cuda_manager.setup(numWorkGroups);
        cuda_manager.execute( cuda_kernel_max<type>, input->get_buffer()->get(), output->get_buffer()->get(), input->size() );
    }
    else
    {
        cuda_manager.setup(numWorkGroups);
        cuda_manager.execute( cuda_kernel_min<type>, input->get_buffer()->get(), output->get_buffer()->get(), input->size() );
    };
    return output;
};

template <typename type>
type vector_cuda<type>::sum()
{
    // Pointer to work
    pointer input = this->copy();
    pointer output = input;

    // Working groups
    auto workGroupSize = cuda_manager.get_threads();
    // auto numWorkGroups = input->size();
    auto numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    // std::cout << "\nwork group size =" << workGroupSize << std::endl;
    // std::cout << "num groups =" << numWorkGroups << std::endl;

    // Main loop
    while(numWorkGroups > workGroupSize)
    {
        output = sum_loop(input, workGroupSize, numWorkGroups);
        input = output;
        // numWorkGroups = input->size();
        numWorkGroups = 1 + ((input->size() - 1)/workGroupSize);   // round up num work groups
    };

    // Finish the reduction work in CPU
    std::vector<type> out = output->std_vector();
    type sum_ = 0;
    for(int i = 0; i < out.size(); i++) sum_ += out[i];
    return sum_;
};

template <typename type>
typename vector_cuda<type>::pointer vector_cuda<type>::sum_loop(typename vector_cuda<type>::pointer input, int workGroupSize, int numWorkGroups)
{
    auto output = vector_cuda<type>::new_pointer(numWorkGroups);

    // cuda_manager.setup(numWorkGroups);
    cuda_manager.setup(input->get_buffer()->size());
    cuda_manager.execute( cuda_kernel_sum<type>, input->get_buffer()->get(), output->get_buffer()->get(), input->get_buffer()->size() );
    return output;
};

// Vectorial dot product. Verify the same number of elements, then product and reduce
template <typename type>
type vector_cuda<type>::dot(const vector_cuda<type> & input)
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
typename vector_cuda<type>::pointer vector_cuda<type>::normalize(type min, type max)
{
    type minv = this->min();
    type maxv = this->max();

    vector_cuda<type>::pointer output = *this - minv;
    output = *(*output*((max - min)/(maxv - minv))) + min;
    return output;
};

template <typename type> template <typename type_cast>
typename vector_cuda<type_cast>::pointer vector_cuda<type>::cast()
{
    int size = this->size();
    auto output = vector_cuda<type_cast>::new_pointer(size);
    
    cuda_manager.setup(size);
    cuda_manager.execute( cuda_kernel_cast<type,type_cast>, buffer->get(), output->get_buffer()->get(), size );
    return output;
};

template <typename type>
void vector_cuda<type>::pad(typename vector_cuda<type>::pointer input, typename vector_cuda<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post)
{
    if(sz.size() == 2)
    {
        cuda_manager.setup(sz);
        cuda_manager.execute( cuda_kernel_pad_2d<type>, input->get_buffer()->get(), output->get_buffer()->get(),
                              pre[0], pre[1], post[0], post[1], sz[0], sz[1] );
        
    }
    else if (sz.size() == 3)
    {
        cuda_manager.setup(sz);
        cuda_manager.execute( cuda_kernel_pad_3d<type>, input->get_buffer()->get(), output->get_buffer()->get(),
                              pre[0], pre[1], pre[2], post[0], post[1], post[2], sz[0], sz[1], sz[2] );
    }
    else ;
};

template <typename type>
void vector_cuda<type>::unpad(typename vector_cuda<type>::pointer input, typename vector_cuda<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post)
{
    if(sz.size() == 2)
    {
        cuda_manager.setup(sz);
        cuda_manager.execute( cuda_kernel_unpad_2d<type>, input->get_buffer()->get(), output->get_buffer()->get(),
                              pre[0], pre[1], post[0], post[1], sz[0], sz[1] );
        
    }
    else if (sz.size() == 3)
    {
        cuda_manager.setup(sz);
        cuda_manager.execute( cuda_kernel_unpad_3d<type>, input->get_buffer()->get(), output->get_buffer()->get(),
                              pre[0], pre[1], pre[2], post[0], post[1], post[2], sz[0], sz[1], sz[2] );
    }
    else ;
};

template <typename type>
void vector_cuda<type>::affine_2d(typename vector_cuda<type>::pointer xi,
                                 typename vector_cuda<type>::pointer yi,
                                 typename vector_cuda<type>::pointer xo,
                                 typename vector_cuda<type>::pointer yo,
                                 typename vector_cuda<type>::pointer p)
{
    int sz = xo->size();
    cuda_manager.setup(sz); // single dim
    cuda_manager.execute( cuda_kernel_affine_2d<type>, xi->get_buffer()->get(), yi->get_buffer()->get(),
                          xo->get_buffer()->get(), yo->get_buffer()->get(),
                          p->get_buffer()->get(), sz);
};

template <typename type>
void vector_cuda<type>::affine_3d(typename vector_cuda<type>::pointer xi,
                                 typename vector_cuda<type>::pointer yi,
                                 typename vector_cuda<type>::pointer zi,
                                 typename vector_cuda<type>::pointer xo,
                                 typename vector_cuda<type>::pointer yo,
                                 typename vector_cuda<type>::pointer zo,
                                 typename vector_cuda<type>::pointer p)
{
    int sz = xo->size();
    cuda_manager.setup(sz); // single dim
    cuda_manager.execute( cuda_kernel_affine_3d<type>, xi->get_buffer()->get(), yi->get_buffer()->get(), zi->get_buffer()->get(),
                            xo->get_buffer()->get(), yo->get_buffer()->get(), zo->get_buffer()->get(),
                            p->get_buffer()->get(), sz);
};

template <typename type>
void vector_cuda<type>::affine_sod_2d(typename vector_cuda<type>::pointer xo,
                                     typename vector_cuda<type>::pointer yo,
                                     typename vector_cuda<type>::pointer xr,
                                     typename vector_cuda<type>::pointer yr,
                                     std::vector<double> & sod)
{
    auto p = vector_cuda<double>::new_pointer(sod.size());
    p->read_ram(sod.data(),sod.size());

    int sz = xo->size();
    cuda_manager.setup(sz); // single dim
    cuda_manager.execute( cuda_kernel_affine_sod_2d<type>, xo->get_buffer()->get(), yo->get_buffer()->get(), 
                         xr->get_buffer()->get(), yr->get_buffer()->get(),
                         p->get_buffer()->get(), sz);
};

template <typename type>
void vector_cuda<type>::affine_sod_3d(typename vector_cuda<type>::pointer xo,
                                     typename vector_cuda<type>::pointer yo,
                                     typename vector_cuda<type>::pointer zo,
                                     typename vector_cuda<type>::pointer xr,
                                     typename vector_cuda<type>::pointer yr,
                                     typename vector_cuda<type>::pointer zr,
                                     std::vector<double> & sod)
{
    auto p = vector_cuda<double>::new_pointer(sod.size());
    p->read_ram(sod.data(),sod.size());

    int sz = xo->size();
    cuda_manager.setup(sz); // single dim
    cuda_manager.execute( cuda_kernel_affine_sod_3d<type>, xo->get_buffer()->get(), yo->get_buffer()->get(), zo->get_buffer()->get(),
                         xr->get_buffer()->get(), yr->get_buffer()->get(), zr->get_buffer()->get(),
                         p->get_buffer()->get(), sz);
};

template <typename type>
void vector_cuda<type>::dfield_2d(typename vector_cuda<type>::pointer xi, typename vector_cuda<type>::pointer yi,
                                 typename vector_cuda<type>::pointer x, typename vector_cuda<type>::pointer y,
                                 typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo)
{
    int sz = xo->size();
    cuda_manager.setup(sz); // single dim
    cuda_manager.execute( cuda_kernel_dfield_2d<type>, xi->get_buffer()->get(), yi->get_buffer()->get(), 
                         x->get_buffer()->get(), y->get_buffer()->get(), 
                         xo->get_buffer()->get(), yo->get_buffer()->get(), sz);
};

template <typename type>
void vector_cuda<type>::dfield_3d(typename vector_cuda<type>::pointer xi, typename vector_cuda<type>::pointer yi, typename vector_cuda<type>::pointer zi,
                                 typename vector_cuda<type>::pointer x, typename vector_cuda<type>::pointer y, typename vector_cuda<type>::pointer z,
                                 typename vector_cuda<type>::pointer xo, typename vector_cuda<type>::pointer yo, typename vector_cuda<type>::pointer zo)
{
    int sz = xo->size();
    cuda_manager.setup(sz); // single dim
    cuda_manager.execute( cuda_kernel_dfield_3d<type>, xi->get_buffer()->get(), yi->get_buffer()->get(), zi->get_buffer()->get(),
                         x->get_buffer()->get(), y->get_buffer()->get(), z->get_buffer()->get(),
                         xo->get_buffer()->get(), yo->get_buffer()->get(), zo->get_buffer()->get(), sz);
};

template <typename type>
std::vector<typename vector_cuda<type>::pointer> vector_cuda<type>::grid2(int w, int h, std::vector<double> & sod)
{
    int size = w*h;
    auto x = vector_cuda<type>::new_pointer(size);
    auto y = vector_cuda<type>::new_pointer(size);
    
    auto p = vector_cuda<double>::new_pointer(sod.size());
    p->read_ram(sod.data(),sod.size());

    std::vector<int> ss({w,h});
    cuda_manager.setup(ss); // multidimensional
    cuda_manager.execute( cuda_kernel_grid_2d<type>, x->get_buffer()->get(), y->get_buffer()->get(), 
                          p->get_buffer()->get(), w, h);
    
    std::vector<typename vector_cuda<type>::pointer> xy(2);
    xy[0] = x;
    xy[1] = y;
    return xy;
};

template <typename type>
std::vector<typename vector_cuda<type>::pointer> vector_cuda<type>::grid3(int w, int h, int l, std::vector<double> & sod)
{
    int size = w*h*l;
    auto x = vector_cuda<type>::new_pointer(size);
    auto y = vector_cuda<type>::new_pointer(size);
    auto z = vector_cuda<type>::new_pointer(size);

    auto p = vector_cuda<double>::new_pointer(sod.size());
    p->read_ram(sod.data(),sod.size());

    std::vector<int> ss({w,h,l});
    cuda_manager.setup(ss); // multidimensional
    cuda_manager.execute( cuda_kernel_grid_3d<type>, x->get_buffer()->get(), y->get_buffer()->get(), 
                          z->get_buffer()->get(), p->get_buffer()->get(), w, h, l);

    std::vector<typename vector_cuda<type>::pointer> xyz(3);
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    return xyz;
};

template <typename type>
void vector_cuda<type>::nearest2( typename vector_cuda<type>::pointer xo, 
                                typename vector_cuda<type>::pointer yo,
                                typename vector_cuda<type>::pointer imgr,
                                typename vector_cuda<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    cuda_manager.setup(out_size); // multidimensional
    cuda_manager.execute( cuda_kernel_nearest_interpolation_2d<type>, 
                        xo->get_buffer()->get(), yo->get_buffer()->get(),
                        imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                        ref_size[0], ref_size[1],
                        out_size[0], out_size[1]);
};

template <typename type>
void vector_cuda<type>::nearest3( typename vector_cuda<type>::pointer xo, 
                                typename vector_cuda<type>::pointer yo,
                                typename vector_cuda<type>::pointer zo,
                                typename vector_cuda<type>::pointer imgr,
                                typename vector_cuda<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    cuda_manager.setup(out_size); // multidimensional
    cuda_manager.execute( cuda_kernel_nearest_interpolation_3d<type>, 
                        xo->get_buffer()->get(), yo->get_buffer()->get(), zo->get_buffer()->get(),
                        imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                        ref_size[0], ref_size[1], ref_size[2],
                        out_size[0], out_size[1], out_size[2]);
};

template <typename type>
void vector_cuda<type>::linear2( typename vector_cuda<type>::pointer xo, 
                                typename vector_cuda<type>::pointer yo,
                                typename vector_cuda<type>::pointer imgr,
                                typename vector_cuda<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    cuda_manager.setup(out_size); // multidimensional
    cuda_manager.execute( cuda_kernel_linear_interpolation_2d<type>, 
                        xo->get_buffer()->get(), yo->get_buffer()->get(),
                        imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                        ref_size[0], ref_size[1],
                        out_size[0], out_size[1]);
};

template <typename type>
void vector_cuda<type>::linear3( typename vector_cuda<type>::pointer xo, 
                                typename vector_cuda<type>::pointer yo,
                                typename vector_cuda<type>::pointer zo,
                                typename vector_cuda<type>::pointer imgr,
                                typename vector_cuda<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    cuda_manager.setup(out_size); // multidimensional
    cuda_manager.execute( cuda_kernel_linear_interpolation_3d<type>, 
                        xo->get_buffer()->get(), yo->get_buffer()->get(), zo->get_buffer()->get(),
                        imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                        ref_size[0], ref_size[1], ref_size[2],
                        out_size[0], out_size[1], out_size[2]);
};

template <typename type>
void vector_cuda<type>::cubic2( typename vector_cuda<type>::pointer xo, 
                                typename vector_cuda<type>::pointer yo,
                                typename vector_cuda<type>::pointer imgr,
                                typename vector_cuda<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    cuda_manager.setup(out_size); // multidimensional
    cuda_manager.execute( cuda_kernel_cubic_interpolation_2d<type>, 
                        xo->get_buffer()->get(), yo->get_buffer()->get(),
                        imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                        ref_size[0], ref_size[1],
                        out_size[0], out_size[1]);
};

template <typename type>
void vector_cuda<type>::cubic3( typename vector_cuda<type>::pointer xo, 
                                typename vector_cuda<type>::pointer yo,
                                typename vector_cuda<type>::pointer zo,
                                typename vector_cuda<type>::pointer imgr,
                                typename vector_cuda<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    cuda_manager.setup(out_size); // multidimensional
    cuda_manager.execute( cuda_kernel_cubic_interpolation_3d<type>, 
                        xo->get_buffer()->get(), yo->get_buffer()->get(), zo->get_buffer()->get(),
                        imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                        ref_size[0], ref_size[1], ref_size[2],
                        out_size[0], out_size[1], out_size[2]);
};

template <typename type>
void vector_cuda<type>::fft(std::vector<pointer> & input, std::vector<pointer> & output, std::vector<int> sz, bool forward)
{
    // std::cout << "CUDA FFT not implemented yet" << std::endl;
    if(sz.size() == 2)
    {
        cuda_manager.setup(sz);
        cuda_manager.execute( cuda_kernel_fft_2d<type>, input[0]->get_buffer()->get(),  input[1]->get_buffer()->get(),
                            output[0]->get_buffer()->get(), output[1]->get_buffer()->get(), sz[0], sz[1], forward);
    }
    else if (sz.size() == 3)
    {
        cuda_manager.setup(sz);
        cuda_manager.execute( cuda_kernel_fft_3d<type>, input[0]->get_buffer()->get(),  input[1]->get_buffer()->get(),
                            output[0]->get_buffer()->get(), output[1]->get_buffer()->get(), sz[0], sz[1], sz[2], forward);
    }
    else ;
};

template <typename type>
void vector_cuda<type>::gradientx( typename vector_cuda<type>::pointer imgr,
                                  typename vector_cuda<type>::pointer imgo,
                                  std::vector<int> ref_size)
{
    cuda_manager.setup(ref_size); // multidimensional
    if (ref_size.size() == 2)
        cuda_manager.execute( cuda_kernel_gradientx<type>, imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                              ref_size[0], ref_size[1], 0 );
    else if (ref_size.size() == 3)
        cuda_manager.execute( cuda_kernel_gradientx<type>, imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                              ref_size[0], ref_size[1], ref_size[2] );
    else ;
};

template <typename type>
void vector_cuda<type>::gradienty( typename vector_cuda<type>::pointer imgr,
                                  typename vector_cuda<type>::pointer imgo,
                                  std::vector<int> ref_size)
{
    cuda_manager.setup(ref_size); // multidimensional
    if (ref_size.size() == 2)
        cuda_manager.execute( cuda_kernel_gradienty<type>, imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                              ref_size[0], ref_size[1], 0 );
    else if (ref_size.size() == 3)
        cuda_manager.execute( cuda_kernel_gradienty<type>, imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                              ref_size[0], ref_size[1], ref_size[2] );
    else ;
};

template <typename type>
void vector_cuda<type>::gradientz( typename vector_cuda<type>::pointer imgr,
                                  typename vector_cuda<type>::pointer imgo,
                                  std::vector<int> ref_size)
{
    cuda_manager.setup(ref_size); // multidimensional
    if (ref_size.size() == 3)
        cuda_manager.execute( cuda_kernel_gradientz<type>, imgr->get_buffer()->get(), imgo->get_buffer()->get(),
                              ref_size[0], ref_size[1], ref_size[2] );
    else ;
};

template <typename type>
void vector_cuda<type>::convolution( typename vector_cuda<type>::pointer imgr,
                                     typename vector_cuda<type>::pointer kernel,
                                     typename vector_cuda<type>::pointer imgo,
                                     std::vector<int> ref_size, std::vector<int> kernel_size)
{
    if (ref_size.size() == 2)
    {
        cuda_manager.setup(ref_size); // multidimensional
        cuda_manager.execute( cuda_kernel_convolution_2d<type>, imgr->get_buffer()->get(), kernel->get_buffer()->get(),
                            imgo->get_buffer()->get(), ref_size[0], ref_size[1], kernel_size[0], kernel_size[1]);
    }
    else if (ref_size.size() == 3)
    {
        cuda_manager.setup(ref_size); // multidimensional
        cuda_manager.execute( cuda_kernel_convolution_3d<type>, imgr->get_buffer()->get(), kernel->get_buffer()->get(), 
                            imgo->get_buffer()->get(), ref_size[0], ref_size[1], ref_size[2], kernel_size[0], kernel_size[1], kernel_size[2]);
    }
    else ;
};

}; //end namespace

#endif