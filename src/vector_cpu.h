/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __VECTOR_CPU_H__
#define __VECTOR_CPU_H__

// std libs
#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <random>       // std::random
#include <cassert>      // assert
#include <cmath>        // math functions
#include <algorithm>    // std::max std::min
// #include <complex>      // std::complex

// local libs
#include "inherit.h"
#include "object.h"
#include "ram_buffer.h"

// fftw
#ifdef IMART_WITH_FFTW
#include <fftw3.h>      // fft library
#endif

// openmp
#ifdef IMART_WITH_OPENMP
#include <omp.h>
#endif

namespace imart
{

// Class object
template <typename type>
class vector_cpu: public inherit<vector_cpu<type>, ram_buffer<type>>
// class vector_cpu: public inherit<vector_cpu<type>, object>, std::vector<type>
{
public:
    //Type definitions
    using self    = vector_cpu;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using object::class_name;
    using object::get_name;
    using ram_buffer<type>::size;
    using ram_buffer<type>::data;
    using ram_buffer<type>::operator[];

    using inherit<vector_cpu<type>, ram_buffer<type>>::inherit;

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
    vector_cpu(): inherit<vector_cpu<type>, ram_buffer<type>>() { class_name = "vector_cpu"; };          // constructor empty
    vector_cpu(int s): inherit<vector_cpu<type>, ram_buffer<type>>(s) { class_name = "vector_cpu"; };    // constructor
    vector_cpu(int s, type value): inherit<vector_cpu<type>, ram_buffer<type>>(s) { assign(value); class_name = "vector_cpu"; };
    vector_cpu(std::initializer_list<type> list): inherit<vector_cpu<type>, ram_buffer<type>>(list) { class_name = "vector_cpu"; };
    vector_cpu(const vector_cpu & input);       // constructor clone
    ~vector_cpu();                              // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const vector_cpu & input);  // copy everything
    virtual void copy_ (const vector_cpu & input);  // share data
    virtual void mimic_(const vector_cpu & input);  // copy meta data

    // ===========================================
    // Get Functions
    // ===========================================
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
    void assert_size(const vector_cpu<type> & input);

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
    vector_cpu<type> & operator = (const vector_cpu & input);
    vector_cpu<type> & operator = (std::initializer_list<type> il);
    pointer operator = (const vector_cpu<type>::pointer input);
    pointer operator + (const vector_cpu<type> & input);
    pointer operator - (const vector_cpu<type> & input);
    pointer operator * (const vector_cpu<type> & input);
    pointer operator / (const vector_cpu<type> & input);
    pointer operator ^ (const vector_cpu<type> & input);

    // Scalar to vector
    pointer operator + (type scalar);
    pointer operator - (type scalar);
    pointer operator * (type scalar);
    pointer operator / (type scalar);
    pointer operator ^ (type scalar);

    // Friend classes to support scalar to vector left hand side
    template<typename type_>
    friend typename vector_cpu<type_>::pointer operator + (type_ scalar, vector_cpu<type_> & input);
    template<typename type_>
    friend typename vector_cpu<type_>::pointer operator - (type_ scalar, const vector_cpu<type_> & input);
    template<typename type_>
    friend typename vector_cpu<type_>::pointer operator * (type_ scalar, vector_cpu<type_> & input);
    template<typename type_>
    friend typename vector_cpu<type_>::pointer operator / (type_ scalar, const vector_cpu<type_> & input);
    template<typename type_>
    friend typename vector_cpu<type_>::pointer operator ^ (type_ scalar, const vector_cpu<type_> & input);

    // Conditions 
    pointer operator == (const vector_cpu<type> & input);
    pointer operator > (const vector_cpu<type> & input);
    pointer operator < (const vector_cpu<type> & input);
    pointer operator >= (const vector_cpu<type> & input);
    pointer operator <= (const vector_cpu<type> & input);

    pointer operator == (type scalar);
    pointer operator > (type scalar);
    pointer operator < (type scalar); // std::shared_ptr<vector_cpu<type>>
    pointer operator >= (type scalar);
    pointer operator <= (type scalar);

    void replace(const vector_cpu<type> & idxs, const vector_cpu<type> & input);
    void replace(const vector_cpu<type> & idxs, type value);
    
    // ===========================================
    // Reduction functions
    // ===========================================
    type min();
    type max();
    type sum();
    // type prod();    // may produce overflow error
    type dot(const vector_cpu<type> & input);

    // ===========================================
    // Functions
    // ===========================================
    pointer normalize(type min = 0.0, type max = 1.0);

    template<typename type_cast>
    typename vector_cpu<type_cast>::pointer cast();

    static void pad(typename vector_cpu<type>::pointer input, typename vector_cpu<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post);
    static void unpad(typename vector_cpu<type>::pointer input, typename vector_cpu<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post);

    static void affine_2d(typename vector_cpu<type>::pointer xi, typename vector_cpu<type>::pointer yi,
                          typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo,
                          typename vector_cpu<type>::pointer p);

    static void affine_3d(typename vector_cpu<type>::pointer xi, typename vector_cpu<type>::pointer yi, typename vector_cpu<type>::pointer zi,
                          typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo, typename vector_cpu<type>::pointer zo,
                          typename vector_cpu<type>::pointer p);

    static void affine_sod_2d(typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo,
                              typename vector_cpu<type>::pointer xr, typename vector_cpu<type>::pointer yr,
                              std::vector<double> & sod);
    static void affine_sod_3d(typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo, typename vector_cpu<type>::pointer zo,
                              typename vector_cpu<type>::pointer xr, typename vector_cpu<type>::pointer yr, typename vector_cpu<type>::pointer zr,
                              std::vector<double> & sod);

    static void  dfield_2d(typename vector_cpu<type>::pointer xi, typename vector_cpu<type>::pointer yi,
                           typename vector_cpu<type>::pointer x, typename vector_cpu<type>::pointer y,
                           typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo);

    static void  dfield_3d(typename vector_cpu<type>::pointer xi, typename vector_cpu<type>::pointer yi, typename vector_cpu<type>::pointer zi,
                           typename vector_cpu<type>::pointer x, typename vector_cpu<type>::pointer y, typename vector_cpu<type>::pointer z,
                           typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo, typename vector_cpu<type>::pointer zo);

    static std::vector<typename vector_cpu<type>::pointer> grid2(int w, int h, std::vector<double> & sod);
    static std::vector<typename vector_cpu<type>::pointer> grid3(int w, int h, int l, std::vector<double> & sod);

    static void nearest2( typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo,
                   typename vector_cpu<type>::pointer imgr, typename vector_cpu<type>::pointer imgo,
                   std::vector<int> ref_size, std::vector<int> out_size);
    static void nearest3( typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo, typename vector_cpu<type>::pointer zo,
                   typename vector_cpu<type>::pointer imgr, typename vector_cpu<type>::pointer imgo,
                   std::vector<int> ref_size, std::vector<int> out_size);

    static void linear2( typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo,
                   typename vector_cpu<type>::pointer imgr, typename vector_cpu<type>::pointer imgo,
                   std::vector<int> ref_size, std::vector<int> out_size);
    static void linear3( typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo, typename vector_cpu<type>::pointer zo,
                   typename vector_cpu<type>::pointer imgr, typename vector_cpu<type>::pointer imgo,
                   std::vector<int> ref_size, std::vector<int> out_size);

    static void linear2_test( typename vector_cpu<type>::pointer xo, 
                                typename vector_cpu<type>::pointer yo,
                                typename vector_cpu<type>::pointer imgr,
                                typename vector_cpu<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<double> & sod, type defaultv);

    static void cubic2( typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo,
                        typename vector_cpu<type>::pointer imgr, typename vector_cpu<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);
    static void cubic3( typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo, typename vector_cpu<type>::pointer zo,
                        typename vector_cpu<type>::pointer imgr, typename vector_cpu<type>::pointer imgo,
                        std::vector<int> ref_size, std::vector<int> out_size);

    #ifdef IMART_WITH_FFTW
    static void fft(std::vector<pointer> & input, std::vector<pointer> & output, std::vector<int> size, bool forward);
    #endif

    static void gradientx( typename vector_cpu<type>::pointer imgr,
                           typename vector_cpu<type>::pointer imgo,
                           std::vector<int> ref_size);

    static void gradienty( typename vector_cpu<type>::pointer imgr,
                           typename vector_cpu<type>::pointer imgo,
                           std::vector<int> ref_size);

    static void gradientz( typename vector_cpu<type>::pointer imgr,
                           typename vector_cpu<type>::pointer imgo,
                           std::vector<int> ref_size);

    static void convolution( typename vector_cpu<type>::pointer imgr,
                             typename vector_cpu<type>::pointer kernel,
                             typename vector_cpu<type>::pointer imgo,
                             std::vector<int> ref_size, std::vector<int> kernel_size);

};

// ===========================================
//          Functions of Class vector_cpu
// ===========================================

// ===========================================
// Constructors
// ===========================================
template <typename type>
vector_cpu<type>::vector_cpu(const vector_cpu<type> & input)
{
    class_name = "vector_cpu";
    clone_(input);          // call the virtual
};

//! Destructor
template <typename type>
vector_cpu<type>::~vector_cpu()
{
    ;
};

// ===========================================
// Create Functions
// ===========================================
template <typename type>
void vector_cpu<type>::clone_(const vector_cpu<type> & input)
{
    // std::cout << "clone_";
    mimic_(input);
    int size = input.size();

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p1[k] = p2[k];
    };
};

template <typename type>
void vector_cpu<type>::copy_(const vector_cpu<type> & input)
{
    // std::cout << "copy_";
    mimic_(input);
};

template <typename type>
void vector_cpu<type>::mimic_(const vector_cpu<type> & input)
{
    // std::cout << "mimic_";
    int size = input.size();
    this->clear();
    this->resize(size);
};

// ===========================================
// Get Functions
// ===========================================
template <typename type>
std::vector<type> vector_cpu<type>::std_vector()
{
    std::vector<type> output(this->size());
    write_ram(output.data(), this->size());
    return output;
    // return static_cast<std::vector<type>>(*this);
};

// ===========================================
// Memory Functions
// ===========================================
template <typename type>
void vector_cpu<type>::read_ram(type * p, int s, int offset)
{
    type * d = this->data() + offset;
    #pragma omp parallel for if (s > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<s; k++)
    {
        *(d+k) = *(p+k);
    };
};

template <typename type>
void vector_cpu<type>::write_ram(type * p, int s, int offset)
{
    type * d = this->data() + offset;
    #pragma omp parallel for if (s > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<s; k++)
    {
        *(p+k) = *(d+k);
    };
};

template <typename type>
void vector_cpu<type>::equal(pointer input)
{
    assert_size(*input);
    type * d = this->data();
    type * p = input->data();
    int n = input->size();
    #pragma omp parallel for if (n > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<n; k++)
    {
        *(d+k) = *(p+k);
    };
};

template <typename type>
void vector_cpu<type>::to_cpu()
{
    ; // do nothing, already on cpu
};

// ===========================================
// Print Functions
// ===========================================
template <typename type>
std::string vector_cpu<type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Vector CPU Information";
    if (msg != "") { title = msg; };

    // Summary of the object information
    ss << object::info(title);
    ss << "Size: \t\t\t" << this->size() << std::endl;
    // ss << "Capacity: \t\t" << this->capacity() << std::endl;
    return ss.str();
};

template <typename type>
std::string vector_cpu<type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; };

    type * p = this->data();
    ss << "[ ";
    for(int k=0; k<this->size(); k++) { ss << p[k] << " "; };
    ss << "]" << std::endl;
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type>
void vector_cpu<type>::assert_size(const vector_cpu<type> & input)
{
    imart_assert(this->size() == input.size(), "Mismatch of vector size");
};

// ===========================================
// Initialization Functions
// ===========================================
template <typename type>
void vector_cpu<type>::zeros()
{
    assign( (type)0.0 );
    // int size = this->size();
    // type * p = this->data();
    // for(int k=0; k<size; k++)
    // {
    //     p[k] = (type)0; // casting to pixel_type
    // };
};

template <typename type>
void vector_cpu<type>::ones()
{
    assign( (type)1.0 );
    // int size = this->size();
    // type * p = this->data();
    // for(int k=0; k<size; k++)
    // {
    //     p[k] = (type)1.0; // casting to pixel_type
    // };
};

// Set all pixel to a fixed value
template <typename type>
void vector_cpu<type>::assign(type value)
{
    // Create pointer
    int size = this->size();
    type * p = this->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p[k] = value;
    };
};

template <typename type>
void vector_cpu<type>::random(float min, float max)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    // std::default_random_engine gen(rd()); //Standard random generator()
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> uniform(min, max);
    int size = this->size();
    type * p = this->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p[k] = (type)uniform(gen); // casting to pixel_type
    };
};

// ===========================================
// Overloading operators
// ===========================================
template <typename type>
vector_cpu<type> & vector_cpu<type>::operator = (const vector_cpu & input)
{
    copy_(input);
    return *this;
};

template <typename type>
vector_cpu<type> & vector_cpu<type>::operator = (std::initializer_list<type> list)
{
    this->clear();
    this->resize(list.size());

    type * d = this->data();
    const type * p = list.begin();
    for(int k=0; k<list.size(); k++)
    {
        *(d+k) = *(p+k);
    };
    return *this;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator = (const vector_cpu<type>::pointer input)
{
    return input;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator + (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic();

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = p1[k] + p2[k];
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator - (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic();

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = p1[k] - p2[k];
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator * (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic();

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = p1[k] * p2[k];
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator / (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic();

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = p1[k] / p2[k];
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator ^ (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic();

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = pow(p1[k],p2[k]);
    };
    return output;
};

// Scalar right hand side
template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator + (type scalar)
{
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic(); // init a image with same poperties
    type * p1 = this->data();
    type * p2 = output->data();
    
    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = p1[k] + scalar;
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator - (type scalar)
{
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic(); // init a image with same poperties
    type * p1 = this->data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = p1[k] - scalar;
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator * (type scalar)
{
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic(); // init a image with same poperties
    type * p1 = this->data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = p1[k] * scalar;
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator / (type scalar)
{
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic(); // init a image with same poperties
    type * p1 = this->data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = p1[k] / scalar;
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator ^ (type scalar)
{
    int size = this->size();
    vector_cpu<type>::pointer output = this->mimic(); // init a image with same poperties
    type * p1 = this->data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = pow(p1[k],scalar);
    };
    return output;
};

// Scalar left hand side
template <typename type>
typename vector_cpu<type>::pointer operator + (type scalar, vector_cpu<type> & input)
{
    return input + scalar;
};

template <typename type>
typename vector_cpu<type>::pointer operator - (type scalar, const vector_cpu<type> & input)
{
    int size = input.size();
    typename vector_cpu<type>::pointer output = input.mimic(); // init a image with same poperties
    type * p1 = input.data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = scalar - p1[k];
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer operator * (type scalar, vector_cpu<type> & input)
{
    return input * scalar;
};

template <typename type>
typename vector_cpu<type>::pointer operator / (type scalar, const vector_cpu<type> & input)
{
    int size = input.size();
    typename vector_cpu<type>::pointer output = input.mimic(); // init a image with same poperties
    type * p1 = input.data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = scalar / p1[k];
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer operator ^ (type scalar, const vector_cpu<type> & input)
{
    int size = input.size();
    typename vector_cpu<type>::pointer output = input.mimic(); // init a image with same poperties
    type * p1 = input.data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = pow(scalar, p1[k]);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator == (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size);

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = (p1[k] == p2[k]);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator > (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size);

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = (p1[k] > p2[k]);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator < (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size);

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = (p1[k] < p2[k]);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator >= (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size);

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = (p1[k] >= p2[k]);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator <= (const vector_cpu<type> & input)
{
    assert_size(input);
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size);

    // Create pointers
    type * p1 = this->data();
    type * p2 = input.data();
    type * p3 = output->data();

    #pragma omp parallel for if (size >= IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p3[k] = (p1[k] <= p2[k]);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator == (type scalar)
{
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size); // init a image with same size
    type * p1 = this->data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = (p1[k] == scalar);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator > (type scalar)
{
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size); // init a image with same size
    type * p1 = this->data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = (p1[k] > scalar);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator < (type scalar)
{
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size); // init a image with same size
    type * p1 = this->data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = (p1[k] < scalar);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator >= (type scalar)
{
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size); // init a image with same size
    type * p1 = this->data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = (p1[k] >= scalar);
    };
    return output;
};

template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::operator <= (type scalar)
{
    int size = this->size();
    auto output = vector_cpu<type>::new_pointer(size); // init a image with same size
    type * p1 = this->data();
    type * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = (p1[k] <= scalar);
    };
    return output;
};

template <typename type>
void vector_cpu<type>::replace(const vector_cpu<type> & idxs, const vector_cpu<type> & input)
{
    assert_size(idxs);
    assert_size(input);
    int size = this->size();

    // Create pointers
    type * p1 = idxs.data();
    type * p2 = this->data();
    type * p3 = input->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        if(p1[k]) p2[k] = p3[k];
    };
    return;
};

template <typename type>
void vector_cpu<type>::replace(const vector_cpu<type> & idxs, type value)
{
    assert_size(idxs);
    int size = this->size();

    // Create pointers
    type * p1 = idxs.data();
    type * p2 = this->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        if(p1[k]) p2[k] = value;
    };
    return;
};

// ===========================================
// Reduction Functions
// ===========================================
template <typename type>
type vector_cpu<type>::min()
{
    type x = 0;
    type * p1 = this->data();
    int size = this->size();

    if (size > 0) { x = p1[0]; };

    #pragma omp parallel for reduction(min:x)
    for(int k=1; k<size; k++)
    {
        x = std::min(x,p1[k]);
    };
    return x;
};

template <typename type>
type vector_cpu<type>::max()
{
    type x = 0;
    type * p1 = this->data();
    int size = this->size();

    if (size > 0) { x = p1[0]; };

    #pragma omp parallel for reduction(max:x)
    for(int k=1; k<size; k++)
    {
        x = std::max(x,p1[k]);
    };
    return x;
};

template <typename type>
type vector_cpu<type>::sum()
{
    type x = 0;
    type * p1 = this->data();
    int size = this->size();    

    #pragma omp parallel for reduction(+:x)
    for(int k=0; k<size; k++)
    {
        x += p1[k];
    };
    return x;
};

// template <typename type>
// type vector_cpu<type>::prod()
// {
//     type x = 1;
//     type * p1 = this->data();
//     int size = this->size();

//     for(int k=0; k<size; k++)
//     {
//         x *= p1[k];
//     };
//     return x;
// };

// Vectorial dot product. Verify the same number of elements, then product and reduce
template <typename type>
type vector_cpu<type>::dot(const vector_cpu<type> & input)
{
    assert_size(input);

    type x = 0;
    type * p1 = this->data();
    type * p2 = input.data();
    int size = this->size();

    #pragma omp parallel for reduction(+:x)
    for(int k=0; k<size; k++)
    {
        x += p1[k]*p2[k];
    };
    return x;
};

// ===========================================
// Functions
// ===========================================
template <typename type>
typename vector_cpu<type>::pointer vector_cpu<type>::normalize(type min, type max)
{
    type minv = this->min();
    type maxv = this->max();

    vector_cpu<type>::pointer output = *this - minv;
    output = *(*output*((max - min)/(maxv - minv))) + min;
    return output;
};

template <typename type> template <typename type_cast> 
typename vector_cpu<type_cast>::pointer vector_cpu<type>::cast()
{
    int size = this->size();
    auto output = vector_cpu<type_cast>::new_pointer(size);
    
    type * p1 = this->data();
    type_cast * p2 = output->data();

    #pragma omp parallel for if (size > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k=0; k<size; k++)
    {
        p2[k] = (type_cast)p1[k];
    };
    return output;
};

template <typename type>
void vector_cpu<type>::pad(typename vector_cpu<type>::pointer input, typename vector_cpu<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post)
{
    type * p1 = input->data();
    type * p2 = output->data();
    if (sz.size() == 2)
    {
        int w = sz[0]; int h = sz[1];
        int start0 = pre[0]; int start1 = pre[1];
        int end0 = post[0]; int end1 = post[1];
        int wo = w + start0 + end0;

        #pragma omp parallel for collapse(2) if (w*h > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int j = 0; j < h; j++)
        {
            for(int i = 0; i < w; i++)
            {
                p2[start0+i+(start1+j)*wo] = p1[i + j*w];
            };
        };

    }
    else if (sz.size() == 3)
    {
        int w = sz[0]; int h = sz[1]; int l = sz[2];
        int start0 = pre[0]; int start1 = pre[1]; int start2 = pre[2];
        int end0 = post[0]; int end1 = post[1]; int end2 = post[2];
        int wo = w + start0 + end0;
        int ho = h + start1 + end1;

        #pragma omp parallel for collapse(3) if (w*h*l > IMART_OPENMP_VECTOR_MIN_SIZE)
        for (int k = 0; k < l; k++)
        {
            for(int j = 0; j < h; j++)
            {
                for(int i = 0; i < w; i++)
                {
                    p2[start0+i+(start1+j)*wo+(start2+k)*wo*ho] = p1[i + j*w + k*w*h];
                };
            };
        };

    }
    else ;
};

template <typename type>
void vector_cpu<type>::unpad(typename vector_cpu<type>::pointer input, typename vector_cpu<type>::pointer output, std::vector<int> sz, std::vector<int> & pre, std::vector<int> & post)
{
    type * p1 = input->data();
    type * p2 = output->data();
    int c = 0;
    if (sz.size() == 2)
    {
        int w = sz[0]; int h = sz[1];
        int start0 = pre[0]; int start1 = pre[1];
        int end0 = post[0]; int end1 = post[1];
        int wo = w + start0 + end0;

        #pragma omp parallel for collapse(2) if (w*h > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int j = 0; j < h; j++)
        {
            for(int i = 0; i < w; i++)
            {
                p2[i + j*w] = p1[start0+i+(start1+j)*wo];
            };
        };
    }
    else if (sz.size() == 3)
    {
        int w = sz[0]; int h = sz[1]; int l = sz[2];
        int start0 = pre[0]; int start1 = pre[1]; int start2 = pre[2];
        int end0 = post[0]; int end1 = post[1]; int end2 = post[2];
        int wo = w + start0 + end0;
        int ho = h + start1 + end1;

        #pragma omp parallel for collapse(3) if (w*h*l > IMART_OPENMP_VECTOR_MIN_SIZE)
        for (int k = 0; k < l; k++)
        {
            for(int j = 0; j < h; j++)
            {
                for(int i = 0; i < w; i++)
                {
                    p2[i + j*w + k*w*h] = p1[start0+i+(start1+j)*wo+(start2+k)*wo*ho];
                };
            };
        };
    }
    else ;
};

template <typename type>
void vector_cpu<type>::affine_2d(typename vector_cpu<type>::pointer xi,
                                 typename vector_cpu<type>::pointer yi,
                                 typename vector_cpu<type>::pointer xo,
                                 typename vector_cpu<type>::pointer yo,
                                 typename vector_cpu<type>::pointer p)
{
    // raw pointers
    type * px = xi->data();
    type * py = yi->data();
    type * pxo = xo->data();
    type * pyo = yo->data();
    type * pp = p->data();

    type a0 = pp[0]; type a1 = pp[1];
    type a2 = pp[2]; type a3 = pp[3];
    type t0 = pp[4]; type t1 = pp[5];

    int sz = xo->size();

    #pragma omp parallel for if (sz > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int i = 0; i < sz; i++)
    {
        pxo[i] = a0*px[i] + a1*py[i] + t0;
        pyo[i] = a2*px[i] + a3*py[i] + t1;
    }
};

template <typename type>
void vector_cpu<type>::affine_3d(typename vector_cpu<type>::pointer xi,
                                 typename vector_cpu<type>::pointer yi,
                                 typename vector_cpu<type>::pointer zi,
                                 typename vector_cpu<type>::pointer xo,
                                 typename vector_cpu<type>::pointer yo,
                                 typename vector_cpu<type>::pointer zo,
                                 typename vector_cpu<type>::pointer p)
{
    // raw pointers
    type * px = xi->data();
    type * py = yi->data();
    type * pz = zi->data();

    type * pxo = xo->data();
    type * pyo = yo->data();
    type * pzo = zo->data();

    type * pp = p->data();
    type a0 = pp[0]; type a1 = pp[1]; type a2 = pp[2];
    type a3 = pp[3]; type a4 = pp[4]; type a5 = pp[5];
    type a6 = pp[6]; type a7 = pp[7]; type a8 = pp[8];
    type t0 = pp[9]; type t1 = pp[10]; type t2 = pp[11];

    int sz = xo->size();

    #pragma omp parallel for if (sz > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int i = 0; i < sz; i++)
    {
        pxo[i] = a0*px[i] + a1*py[i] + a2*pz[i] + t0;
        pyo[i] = a3*px[i] + a4*py[i] + a5*pz[i] + t1;
        pzo[i] = a6*px[i] + a7*py[i] + a8*pz[i] + t2;
    };
};

template <typename type>
void vector_cpu<type>::affine_sod_2d(typename vector_cpu<type>::pointer xo,
                                     typename vector_cpu<type>::pointer yo,
                                     typename vector_cpu<type>::pointer xr,
                                     typename vector_cpu<type>::pointer yr,
                                     std::vector<double> & sod)
{
    // raw pointers
    type * px = xo->data();
    type * py = yo->data();
    type * pxo = xr->data();
    type * pyo = yr->data();

    double s0 = sod[0]; double s1 = sod[1];
    double o0 = sod[2]; double o1 = sod[3];
    double d0 = sod[4]; double d1 = sod[5];
    double d2 = sod[6]; double d3 = sod[7];

    int sz = xo->size();

    #pragma omp parallel for if (sz > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int i = 0; i < sz; i++)
    {
        pxo[i] = d0*s0*px[i] + d1*s1*py[i] + o0;
        pyo[i] = d2*s0*px[i] + d3*s1*py[i] + o1;
    }
};

template <typename type>
void vector_cpu<type>::affine_sod_3d(typename vector_cpu<type>::pointer xo,
                                     typename vector_cpu<type>::pointer yo,
                                     typename vector_cpu<type>::pointer zo,
                                     typename vector_cpu<type>::pointer xr,
                                     typename vector_cpu<type>::pointer yr,
                                     typename vector_cpu<type>::pointer zr,
                                     std::vector<double> & sod)
{
    // raw pointers
    type * px = xo->data();
    type * py = yo->data();
    type * pz = zo->data();

    type * pxo = xr->data();
    type * pyo = yr->data();
    type * pzo = zr->data();

    double s0 = sod[0]; double s1 = sod[1]; double s2 = sod[2];
    double o0 = sod[3]; double o1 = sod[4]; double o2 = sod[5];
    double d0 = sod[6]; double d1 = sod[7]; double d2 = sod[8];
    double d3 = sod[9]; double d4 = sod[10]; double d5 = sod[11];
    double d6 = sod[12]; double d7 = sod[13]; double d8 = sod[14];

    int sz = xo->size();

    #pragma omp parallel for if (sz > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int i = 0; i < sz; i++)
    {
        pxo[i] = d0*s0*px[i] + d1*s1*py[i] + d2*s2*pz[i] + o0;
        pyo[i] = d3*s0*px[i] + d4*s1*py[i] + d5*s2*pz[i] + o1;
        pzo[i] = d6*s0*px[i] + d7*s1*py[i] + d8*s2*pz[i] + o2;
    };
};

template <typename type>
void vector_cpu<type>::dfield_2d(typename vector_cpu<type>::pointer xi, typename vector_cpu<type>::pointer yi,
                                 typename vector_cpu<type>::pointer x, typename vector_cpu<type>::pointer y, 
                                 typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo)
{
    // raw pointers
    type * pxi = xi->data();
    type * pyi = yi->data();

    type * px = x->data();
    type * py = y->data();

    type * pxo = xo->data();
    type * pyo = yo->data();

    int sz = xo->size();

    #pragma omp parallel for if (sz > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int i = 0; i < sz; i++)
    {
        pxo[i] = pxi[i] + px[i];
        pyo[i] = pyi[i] + py[i];
    }
};

template <typename type>
void vector_cpu<type>::dfield_3d(typename vector_cpu<type>::pointer xi, typename vector_cpu<type>::pointer yi, typename vector_cpu<type>::pointer zi,
                                 typename vector_cpu<type>::pointer x, typename vector_cpu<type>::pointer y, typename vector_cpu<type>::pointer z,
                                 typename vector_cpu<type>::pointer xo, typename vector_cpu<type>::pointer yo, typename vector_cpu<type>::pointer zo)
{
    // raw pointers
    type * pxi = xi->data();
    type * pyi = yi->data();
    type * pzi = zi->data();

    type * px = x->data();
    type * py = y->data();
    type * pz = z->data();

    type * pxo = xo->data();
    type * pyo = yo->data();
    type * pzo = zo->data();

    int sz = xi->size();

    #pragma omp parallel for if (sz > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int i = 0; i < sz; i++)
    {
        pxo[i] = pxi[i] + px[i];
        pyo[i] = pyi[i] + py[i];
        pzo[i] = pzi[i] + pz[i];
    };
};

template <typename type>
std::vector<typename vector_cpu<type>::pointer> vector_cpu<type>::grid2(int w, int h, std::vector<double> & sod)
{
    int size = w*h;
    auto x = vector_cpu<type>::new_pointer(size);
    auto y = vector_cpu<type>::new_pointer(size);

    type * px = x->data();
    type * py = y->data();

    double s0 = sod[0]; double s1 = sod[1];
    double o0 = sod[2]; double o1 = sod[3];
    double d0 = sod[4]; double d1 = sod[5];
    double d2 = sod[6]; double d3 = sod[7];

    #pragma omp parallel for collapse(2) if (w*h > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int j = 0; j < h; j++)
    {
        for(int i = 0; i < w; i++)
        {
            // without direction
            // px[i+j*w] = s0*i + o0;
            // py[i+j*w] = s1*j + o1;
            // with direction
            px[i+j*w] = d0*s0*i + d1*s1*j + o0;
            py[i+j*w] = d2*s0*i + d3*s1*j + o1;
        }
    }
       
    std::vector<typename vector_cpu<type>::pointer> xy(2);
    xy[0] = x;
    xy[1] = y;
    return xy;
};

template <typename type>
std::vector<typename vector_cpu<type>::pointer> vector_cpu<type>::grid3(int w, int h, int l, std::vector<double> & sod)
{
    int size = w*h*l;
    auto x = vector_cpu<type>::new_pointer(size);
    auto y = vector_cpu<type>::new_pointer(size);
    auto z = vector_cpu<type>::new_pointer(size);

    type * px = x->data();
    type * py = y->data();
    type * pz = z->data();

    double s0 = sod[0]; double s1 = sod[1]; double s2 = sod[2];
    double o0 = sod[3]; double o1 = sod[4]; double o2 = sod[5];
    double d0 = sod[6]; double d1 = sod[7]; double d2 = sod[8];
    double d3 = sod[9]; double d4 = sod[10]; double d5 = sod[11];
    double d6 = sod[12]; double d7 = sod[13]; double d8 = sod[14];

    #pragma omp parallel for collapse(3) if (w*h*l > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k = 0; k < l; k++)
    {
        for(int j = 0; j < h; j++)
        {
            for(int i = 0; i < w; i++)
            {
                // without direction
                // px[i + j*w + k*w*h] = s0*i + o0;
                // py[i + j*w + k*w*h] = s1*j + o1;
                // pz[i + j*w + k*w*h] = s2*k + o2;
                // with direction
                px[i + j*w + k*w*h] = d0*s0*i + d1*s1*j + d2*s2*k + o0;
                py[i + j*w + k*w*h] = d3*s0*i + d4*s1*j + d5*s2*k + o1;
                pz[i + j*w + k*w*h] = d6*s0*i + d7*s1*j + d8*s2*k + o2;
            };
        };
    };
       
    std::vector<typename vector_cpu<type>::pointer> xyz(3);
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
    return xyz;
};

template <typename type>
void vector_cpu<type>::nearest2( typename vector_cpu<type>::pointer xo, 
                                typename vector_cpu<type>::pointer yo,
                                typename vector_cpu<type>::pointer imgr,
                                typename vector_cpu<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    int n0 = out_size[0]; int n1 = out_size[1];
    int w = ref_size[0]; int h = ref_size[1];

    type * pxo = xo->data();
    type * pyo = yo->data();
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    #pragma omp parallel for collapse(2) if (n0*n1 > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int j = 0; j < n1; j++)
    {
        for(int i = 0; i < n0; i++)
        {
            int x = round(pxo[i + j*n0]);
            int y = round(pyo[i + j*n0]);
            if(x >= 0 && x < w && y >= 0 && y < h)
            {
                pimgo[i + j*n0] = pimgr[x + y*w];
            };
        };
    };
};

template <typename type>
void vector_cpu<type>::nearest3( typename vector_cpu<type>::pointer xo, 
                                typename vector_cpu<type>::pointer yo,
                                typename vector_cpu<type>::pointer zo,
                                typename vector_cpu<type>::pointer imgr,
                                typename vector_cpu<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    int n0 = out_size[0]; int n1 = out_size[1]; int n2 = out_size[1];
    int w = ref_size[0]; int h = ref_size[1]; int l = ref_size[1];

    type * pxo = xo->data();
    type * pyo = yo->data();
    type * pzo = zo->data();
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    #pragma omp parallel for collapse(3) if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k = 0; k < n2; k++)
    {
        for(int j = 0; j < n1; j++)
        {
            for(int i = 0; i < n0; i++)
            {
                int x = round(pxo[i + j*n0 + k*n0*n1]);
                int y = round(pyo[i + j*n0 + k*n0*n1]);
                int z = round(pzo[i + j*n0 + k*n0*n1]);
                if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l)
                {
                    pimgo[i + j*n0 + k*n0*n1] = pimgr[x + y*w + z*w*h];
                };
            };
        };
    };
};

// function
// template <typename type>
// type ilinear2(int idx, type * pxo, type * pyo, type * pimgr, type * pimgo)
// {

// };

// class
template <typename type>
class ilinear2
{
public:
    // variables
    type zero = 0.01;
    int w; int h;
    type * px; type * py; type * ref; type * out;

    ilinear2(int ww, int hh, type * pxo, type * pyo, type * pimgr, type * pimgo)
    {
        w = ww; h = hh;
        px = pxo; py = pyo; ref = pimgr; out = pimgo;
    }

    void interpolate(int idx)
    {
        // type xt = px[idx];
        // int x = floor(xt);
        // if (x < 0 || x > w) return; // outside x

        // type yt = py[idx];
        // int y = floor(yt);
        // if (y < 0 || y > h) return; // outside y

        type xt = px[idx];
        int x = floor(xt);
        type yt = py[idx];
        int y = floor(yt);

        if(x >= 0 && x < w && y >= 0 && y < h)
        {
            type dx = xt - (type)x;
            type dy = yt - (type)y;
            if (dx < zero && dy < zero)
            {
                out[idx] = ref[x+y*w];
            }
            else if (dy < zero || y >= h - 1) // same y
            {
                out[idx] = ref[x+y*w]*(1-dx) + ref[x+1+y*w]*(dx);
            }
            else if (dx < zero || x >= w - 1) // same x
            {
                out[idx] = ref[x+y*w]*(1-dy) + ref[x+(y+1)*w]*(dy);
            }
            else
            {
                // compute case x & y
                type dxdy = dx*dy;
                type r = ref[x+y*w]*(1-dx-dy+dxdy) + ref[x+1+y*w]*(dx-dxdy) + ref[x+(y+1)*w]*(dy-dxdy) + ref[x+1+(y+1)*w]*dxdy;
                out[idx] = r;
            };
        };
    };
};


// optimized
template <typename type>
void vector_cpu<type>::linear2( typename vector_cpu<type>::pointer xo, 
                                typename vector_cpu<type>::pointer yo,
                                typename vector_cpu<type>::pointer imgr,
                                typename vector_cpu<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    int n0 = out_size[0]; int n1 = out_size[1];
    int w = ref_size[0]; int h = ref_size[1];

    type * pxo = xo->data();
    type * pyo = yo->data();
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    ilinear2<type> lin2(w, h, pxo, pyo, pimgr, pimgo);

    #pragma omp parallel for collapse(2) if (n0*n1 > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int j = 0; j < n1; j++)
    {
        for(int i = 0; i < n0; i++)
        {
            int idx = i + j*n0;
            lin2.interpolate(idx);
        };
    };
};

// template <typename type>
// void vector_cpu<type>::linear2( typename vector_cpu<type>::pointer xo, 
//                                 typename vector_cpu<type>::pointer yo,
//                                 typename vector_cpu<type>::pointer imgr,
//                                 typename vector_cpu<type>::pointer imgo,
//                                 std::vector<int> ref_size, std::vector<int> out_size)
// {
//     int n0 = out_size[0]; int n1 = out_size[1];
//     int w = ref_size[0]; int h = ref_size[1];

//     type * pxo = xo->data();
//     type * pyo = yo->data();
//     type * pimgr = imgr->data();
//     type * pimgo = imgo->data();

//     #pragma omp parallel for collapse(2) if (n0*n1 > IMART_OPENMP_VECTOR_MIN_SIZE)
//     for(int j = 0; j < n1; j++)
//     {
//         for(int i = 0; i < n0; i++)
//         {
//             type xt = pxo[i + j*n0];
//             type yt = pyo[i + j*n0];
//             int x = floor(xt);
//             int y = floor(yt);

//             if(x >= 0 && x < w - 1 && y >= 0 && y < h - 1)
//             {
//                 type dx = xt - (type)x;
//                 type dy = yt - (type)y;
//                 type dxdy = dx*dy;
//                 type r = pimgr[x+y*w]*(1-dx-dy+dxdy) + pimgr[x+1+y*w]*(dx-dxdy) + pimgr[x+(y+1)*w]*(dy-dxdy) + pimgr[x+1+(y+1)*w]*dxdy;
//                 pimgo[i + j*n0] = r;
//             }
//             else if(x >= 0 && x < w - 1 && y == h - 1) // border case
//             {
//                 type dx = xt - (type)x;
//                 type r = pimgr[x+y*w]*(1-dx) + pimgr[x+1+y*w]*(dx);
//                 pimgo[i + j*n0] = r;
//             }
//             else if(x >= 0 && x == w - 1 && y < h - 1) // border case
//             {
//                 type dy = yt - (type)y;
//                 type r = pimgr[x+y*w]*(1-dy) + pimgr[x+(y+1)*w]*(dy);
//                 pimgo[i + j*n0] = r;
//             };
//         };
//     };
// };

// template <typename type>
// void vector_cpu<type>::linear2( typename vector_cpu<type>::pointer xo, 
//                                 typename vector_cpu<type>::pointer yo,
//                                 typename vector_cpu<type>::pointer imgr,
//                                 typename vector_cpu<type>::pointer imgo,
//                                 std::vector<int> ref_size, std::vector<int> out_size)
// {
//     int n0 = out_size[0]; int n1 = out_size[1];
//     int w = ref_size[0]; int h = ref_size[1];

//     type * pxo = xo->data();
//     type * pyo = yo->data();
//     type * pimgr = imgr->data();
//     type * pimgo = imgo->data();

//     #pragma omp parallel for collapse(2) if (n0*n1 > IMART_OPENMP_VECTOR_MIN_SIZE)
//     for(size_t j = 0; j < n1; j++)
//     {
//         for(size_t i = 0; i < n0; i++)
//         {
//             size_t p = i + j*n0;
//             type xt = pxo[p];
//             type yt = pyo[p];
//             size_t x = floor(xt);
//             size_t y = floor(yt);

//             if(x >= 0 && x < w && y >= 0 && y < h - 1)
//             {
//                 type dx = xt - (type)x;
//                 type dy = yt - (type)y;
//                 type dxdy = dx*dy;
//                 type r = pimgr[x+y*w]*(1-dx-dy+dxdy) + pimgr[x+1+y*w]*(dx-dxdy) + pimgr[x+(y+1)*w]*(dy-dxdy) + pimgr[x+1+(y+1)*w]*dxdy;
//                 pimgo[p] = r;
//             }
//             else if(x >= 0 && x < w && y == h - 1) // border case
//             {
//                 type dx = xt - (type)x;
//                 type r = pimgr[x+y*w]*(1-dx) + pimgr[x+1+y*w]*(dx);
//                 pimgo[p] = r;
//             };
//         };
//     };
// };


// template <typename type>
// void vector_cpu<type>::linear2_test( typename vector_cpu<type>::pointer xo, 
//                                 typename vector_cpu<type>::pointer yo,
//                                 typename vector_cpu<type>::pointer imgr,
//                                 typename vector_cpu<type>::pointer imgo,
//                                 std::vector<int> ref_size, std::vector<double> & sod, type defaultv)
// {
//     int w = ref_size[0]; int h = ref_size[1];
//     type * px = xo->data();
//     type * py = yo->data();
//     type * pimgr = imgr->data();
//     type * pimgo = imgo->data();

//     double s0 = sod[0]; double s1 = sod[1];
//     double o0 = sod[2]; double o1 = sod[3];
//     double d0 = sod[4]; double d1 = sod[5];
//     double d2 = sod[6]; double d3 = sod[7];

//     int sz = xo->size();

//     #pragma omp parallel for if (sz > IMART_OPENMP_VECTOR_MIN_SIZE)
//     for(size_t i = 0; i < sz; i++)
//     {
//         type xt = d0*s0*px[i] + d1*s1*py[i] + o0;
//         type yt = d2*s0*px[i] + d3*s1*py[i] + o1;
//         size_t x = floor(xt);
//         size_t y = floor(yt);

//         if(x >= 0 && x < w && y >= 0 && y < h - 1)
//         {
//             type dx = xt - (type)x;
//             type dy = yt - (type)y;
//             type dxdy = dx*dy;
//             type r = pimgr[x+y*w]*(1-dx-dy+dxdy) + pimgr[x+1+y*w]*(dx-dxdy) + pimgr[x+(y+1)*w]*(dy-dxdy) + pimgr[x+1+(y+1)*w]*dxdy;
//             pimgo[i] = r;
//         }
//         else if(x >= 0 && x < w && y == h - 1) // border case
//         {
//             type dx = xt - (type)x;
//             type r = pimgr[x+y*w]*(1-dx) + pimgr[x+1+y*w]*(dx);
//             pimgo[i] = r;
//         }
//         else
//         {
//             pimgo[i] = defaultv;
//         };
//     };
// };

// template <typename type>
// void vector_cpu<type>::linear2_test( typename vector_cpu<type>::pointer xo, 
//                                 typename vector_cpu<type>::pointer yo,
//                                 typename vector_cpu<type>::pointer imgr,
//                                 typename vector_cpu<type>::pointer imgo,
//                                 std::vector<int> ref_size, std::vector<double> & sod, type defaultv)
// {
//     size_t w = ref_size[0]; size_t h = ref_size[1];
//     type * px = xo->data();
//     type * py = yo->data();
//     type * pimgr = imgr->data();
//     type * pimgo = imgo->data();

//     double s0 = sod[0]; double s1 = sod[1];
//     double o0 = sod[2]; double o1 = sod[3];
//     double d0 = sod[4]; double d1 = sod[5];
//     double d2 = sod[6]; double d3 = sod[7];

//     int sz = xo->size();

//     #pragma omp parallel for if (sz > IMART_OPENMP_VECTOR_MIN_SIZE)
//     for(size_t i = 0; i < sz; i++)
//     {
//         type xt = d0*s0*px[i] + d1*s1*py[i] + o0;
//         type yt = d2*s0*px[i] + d3*s1*py[i] + o1;
//         size_t x = floor(xt);
//         size_t y = floor(yt);

//         if(x >= 0 && x < w && y >= 0 && y < h - 1)
//         {
//             type dx = xt - (type)x;
//             type dy = yt - (type)y;
//             type dxdy = dx*dy;
//             size_t addr = x+y*w;
//             // type a[4] = {pimgr[addr], pimgr[addr+1], pimgr[addr+w], pimgr[addr+1+w]};
//             // type r = a[0]*(1-dx-dy+dxdy) + a[1]*(dx-dxdy) + a[2]*(dy-dxdy) + a[3]*dxdy;
//             type r = pimgr[addr]*(1-dx-dy+dxdy) + pimgr[addr+1]*(dx-dxdy) + pimgr[addr+w]*(dy-dxdy) + pimgr[addr+1+w]*dxdy;
//             pimgo[i] = r;
//         }
//         else if(x >= 0 && x < w && y == h - 1) // border case
//         {
//             type dx = xt - (type)x;
//             type r = pimgr[x+y*w]*(1-dx) + pimgr[x+1+y*w]*(dx);
//             pimgo[i] = r;
//         }
//         else
//         {
//             pimgo[i] = defaultv;
//         };
//         // pimgo[i] = defaultv;
//     };
// };

template <typename type>
class ilinear3
{
public:
    // variables
    type zero = 0.01;
    int w; int h; int l;
    type * px; type * py; type * pz; type * ref; type * out;

    ilinear3(int ww, int hh, int ll, type * pxo, type * pyo, type * pzo, type * pimgr, type * pimgo)
    {
        w = ww; h = hh; l = ll;
        px = pxo; py = pyo; pz = pzo; ref = pimgr; out = pimgo;
    }

    void interpolate(int idx)
    {
        type xt = px[idx];
        type yt = py[idx];
        type zt = pz[idx];
        int x = floor(xt);
        int y = floor(yt);
        int z = floor(zt);
        if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l)
        {
            type dx = xt - (type)x;
            type dy = yt - (type)y;
            type dz = zt - (type)z;
            if (dx <= zero && dy <= zero && dz <= zero)
            {
                out[idx] = ref[x+y*w+z*w*h];
            }
            else if (dz <= zero || z >= l - 1) // same z
            {   
                if (dy <= zero || y >= h - 1) // same y
                {
                    out[idx] = ref[x+y*w+z*w*h]*(1-dx) + ref[x+1+y*w+z*w*h]*(dx);
                }
                else if (dx <= zero || x >= w - 1) // same x
                {
                    out[idx] = ref[x+y*w+z*w*h]*(1-dy) + ref[x+(y+1)*w+z*w*h]*(dy);
                }
                else
                {
                    // compute case x & y
                    type dxdy = dx*dy;
                    type r = ref[x+y*w+z*w*h]*(1-dx-dy+dxdy) + ref[x+1+y*w+z*w*h]*(dx-dxdy) + ref[x+(y+1)*w+z*w*h]*(dy-dxdy) + ref[x+1+(y+1)*w+z*w*h]*dxdy;
                    out[idx] = r;
                };
            }
            else if (dy <= zero || y >= h - 1) // same y
            {
                if (dx <= zero || x >= w - 1) // same x
                {
                    out[idx] = ref[x+y*w+z*w*h]*(1-dz) + ref[x+y*w+(z+1)*w*h]*(dz);
                    return;
                }
                else
                {
                    // compute case x & z
                    type dxdz = dx*dz;
                    type r = ref[x+y*w+z*w*h]*(1-dx-dz+dxdz) + ref[x+1+y*w+z*w*h]*(dx-dxdz) + ref[x+y*w+(z+1)*w*h]*(dz-dxdz) + ref[x+1+y*w+(z+1)*w*h]*dxdz;
                    out[idx] = r;
                };
            }
            else if (dx <= zero || x >= w - 1) // same x
            {
                // compute case y & z
                type dydz = dy*dz;
                type r = ref[x+y*w+z*w*h]*(1-dy-dz+dydz) + ref[x+(y+1)*w+z*w*h]*(dy-dydz) + ref[x+y*w+(z+1)*w*h]*(dz-dydz) + ref[x+(y+1)*w+(z+1)*w*h]*dydz;
                out[idx] = r;
            }
            else
            {
                // compute case x & y & z
                type dxdy = dx*dy;
                type rv = ref[x+y*w+z*w*h]*(1-dx-dy+dxdy) + ref[x+1+y*w+z*w*h]*(dx-dxdy) + ref[x+(y+1)*w+z*w*h]*(dy-dxdy) + ref[x+1+(y+1)*w+z*w*h]*dxdy;
                type rw = ref[x+y*w+(z+1)*w*h]*(1-dx-dy+dxdy) + ref[x+1+y*w+(z+1)*w*h]*(dx-dxdy) + ref[x+(y+1)*w+(z+1)*w*h]*(dy-dxdy) + ref[x+1+(y+1)*w+(z+1)*w*h]*dxdy;
                type r = rv*(1-dz) + rw*dz;
                out[idx] = r;
            };
        };
    };
};

template <typename type>
void vector_cpu<type>::linear3( typename vector_cpu<type>::pointer xo, 
                                typename vector_cpu<type>::pointer yo,
                                typename vector_cpu<type>::pointer zo,
                                typename vector_cpu<type>::pointer imgr,
                                typename vector_cpu<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    int n0 = out_size[0]; int n1 = out_size[1]; int n2 = out_size[2];
    int w = ref_size[0]; int h = ref_size[1]; int l = ref_size[2];

    type * pxo = xo->data();
    type * pyo = yo->data();
    type * pzo = zo->data();
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    ilinear3<type> lin3(w, h, l, pxo, pyo, pzo, pimgr, pimgo);

    #pragma omp parallel for collapse(3) if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k = 0; k < n2; k++)
    {
        for(int j = 0; j < n1; j++)
        {
            for(int i = 0; i < n0; i++)
            {
                int idx = i + j*n0 + k*n0*n1;
                lin3.interpolate(idx);
            };
        };
    };
};

// template <typename type>
// void vector_cpu<type>::linear3( typename vector_cpu<type>::pointer xo, 
//                                 typename vector_cpu<type>::pointer yo,
//                                 typename vector_cpu<type>::pointer zo,
//                                 typename vector_cpu<type>::pointer imgr,
//                                 typename vector_cpu<type>::pointer imgo,
//                                 std::vector<int> ref_size, std::vector<int> out_size)
// {
//     int n0 = out_size[0]; int n1 = out_size[1]; int n2 = out_size[2];
//     int w = ref_size[0]; int h = ref_size[1]; int l = ref_size[2];

//     type * pxo = xo->data();
//     type * pyo = yo->data();
//     type * pzo = zo->data();
//     type * pimgr = imgr->data();
//     type * pimgo = imgo->data();

//     #pragma omp parallel for collapse(3) if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
//     for(int k = 0; k < n2; k++)
//     {
//         for(int j = 0; j < n1; j++)
//         {
//             for(int i = 0; i < n0; i++)
//             {
//                 type xt = pxo[i + j*n0 + k*n0*n1];
//                 type yt = pyo[i + j*n0 + k*n0*n1];
//                 type zt = pzo[i + j*n0 + k*n0*n1];
//                 int x = floor(xt);
//                 int y = floor(yt);
//                 int z = floor(zt);
//                 if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l - 1)
//                 {
//                     type dx = xt - (type)x;
//                     type dy = yt - (type)y;
//                     type dxdy = dx*dy;
//                     type rv = pimgr[x+y*w+z*w*h]*(1-dx-dy+dxdy) + pimgr[x+1+y*w+z*w*h]*(dx-dxdy) + pimgr[x+(y+1)*w+z*w*h]*(dy-dxdy) + pimgr[x+1+(y+1)*w+z*w*h]*dxdy;
//                     type rw = pimgr[x+y*w+(z+1)*w*h]*(1-dx-dy+dxdy) + pimgr[x+1+y*w+(z+1)*w*h]*(dx-dxdy) + pimgr[x+(y+1)*w+(z+1)*w*h]*(dy-dxdy) + pimgr[x+1+(y+1)*w+(z+1)*w*h]*dxdy;
//                     type dz = zt - (type)z;
//                     type r = rv*(1-dz) + rw*dz;
//                     pimgo[i + j*n0 + k*n0*n1] = r;
//                 }
//                 else if(x >= 0 && x < w && y >= 0 && y < h && z == l - 1) // border case
//                 {
//                     type dx = xt - (type)x;
//                     type dy = yt - (type)y;
//                     type dxdy = dx*dy;
//                     type rv = pimgr[x+y*w+z*w*h]*(1-dx-dy+dxdy) + pimgr[x+1+y*w+z*w*h]*(dx-dxdy) + pimgr[x+(y+1)*w+z*w*h]*(dy-dxdy) + pimgr[x+1+(y+1)*w+z*w*h]*dxdy;
//                     pimgo[i + j*n0 + k*n0*n1] = rv;
//                 };
//             };
//         };
//     };
// };

template <typename type>
type cubic(type p[4], type x) 
{
    return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
};

template <typename type>
void vector_cpu<type>::cubic2(  typename vector_cpu<type>::pointer xo, 
                                typename vector_cpu<type>::pointer yo,
                                typename vector_cpu<type>::pointer imgr,
                                typename vector_cpu<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    int n0 = out_size[0]; int n1 = out_size[1];
    int w = ref_size[0]; int h = ref_size[1];

    type * pxo = xo->data();
    type * pyo = yo->data();
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    #pragma omp parallel for if (n0*n1 > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int j = 0; j < n1; j++)
    {
        for(int i = 0; i < n0; i++)
        {
            type xt = pxo[i + j*n0];
            type yt = pyo[i + j*n0];
            int x = floor(xt);
            int y = floor(yt);

            if(x >= 1 && x < w - 2 && y >= 1 && y < h - 2)
            {
                type dx = xt - (type)x;
                type dy = yt - (type)y;

                type r0[4] = {pimgr[x-1+(y-1)*w], pimgr[x+(y-1)*w], pimgr[x+1+(y-1)*w], pimgr[x+2+(y-1)*w]};
                type r1[4] = {pimgr[x-1+(y)*w]  , pimgr[x+(y)*w]  , pimgr[x+1+(y)*w]  , pimgr[x+2+(y)*w]};
                type r2[4] = {pimgr[x-1+(y+1)*w], pimgr[x+(y+1)*w], pimgr[x+1+(y+1)*w], pimgr[x+2+(y+1)*w]};
                type r3[4] = {pimgr[x-1+(y+2)*w], pimgr[x+(y+2)*w], pimgr[x+1+(y+2)*w], pimgr[x+2+(y+2)*w]};
                
                type r[4] = {cubic(r0, dx), cubic(r1, dx), cubic(r2, dx), cubic(r3, dx) };
                pimgo[i + j*n0] = cubic(r, dy);
            }
            else if(x >= 0 && x < w && y >= 0 && y < h - 1) //linear interpolation otherwise
            {
                type dx = xt - (type)x;
                type dy = yt - (type)y;
                type dxdy = dx*dy;
                type r = pimgr[x+y*w]*(1-dx-dy+dxdy) + pimgr[x+1+y*w]*(dx-dxdy) + pimgr[x+(y+1)*w]*(dy-dxdy) + pimgr[x+1+(y+1)*w]*dxdy;
                pimgo[i + j*n0] = r;
            }
            else if(x >= 0 && x < w && y == h - 1) // border case
            {
                type dx = xt - (type)x;
                type r = pimgr[x+y*w]*(1-dx) + pimgr[x+1+y*w]*(dx);
                pimgo[i + j*n0] = r;
            };
        };
    };
};

template <typename type>
void vector_cpu<type>::cubic3(  typename vector_cpu<type>::pointer xo, 
                                typename vector_cpu<type>::pointer yo,
                                typename vector_cpu<type>::pointer zo,
                                typename vector_cpu<type>::pointer imgr,
                                typename vector_cpu<type>::pointer imgo,
                                std::vector<int> ref_size, std::vector<int> out_size)
{
    int n0 = out_size[0]; int n1 = out_size[1]; int n2 = out_size[2];
    int w = ref_size[0]; int h = ref_size[1]; int l = ref_size[2];

    type * pxo = xo->data();
    type * pyo = yo->data();
    type * pzo = zo->data();
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    #pragma omp parallel for if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
    for(int k = 0; k < n2; k++)
    {
        for(int j = 0; j < n1; j++)
        {
            for(int i = 0; i < n0; i++)
            {
                type xt = pxo[i + j*n0 + k*n0*n1];
                type yt = pyo[i + j*n0 + k*n0*n1];
                type zt = pzo[i + j*n0 + k*n0*n1];
                int x = floor(xt);
                int y = floor(yt);
                int z = floor(zt);
                if(x >= 1 && x < w - 2 && y >= 1 && y < h - 2 && z >= 1 && z < l - 2)
                {
                    type dx = xt - (type)x;
                    type dy = yt - (type)y;
                    type dz = zt - (type)z;

                    type r00[4] = {pimgr[x-1+(y-1)*w+(z-1)*w*h], pimgr[x+(y-1)*w+(z-1)*w*h], pimgr[x+1+(y-1)*w+(z-1)*w*h], pimgr[x+2+(y-1)*w+(z-1)*w*h]};
                    type r01[4] = {pimgr[x-1+(y)*w+(z-1)*w*h]  , pimgr[x+(y)*w+(z-1)*w*h]  , pimgr[x+1+(y)*w+(z-1)*w*h]  , pimgr[x+2+(y)*w+(z-1)*w*h]};
                    type r02[4] = {pimgr[x-1+(y+1)*w+(z-1)*w*h], pimgr[x+(y+1)*w+(z-1)*w*h], pimgr[x+1+(y+1)*w+(z-1)*w*h], pimgr[x+2+(y+1)*w+(z-1)*w*h]};
                    type r03[4] = {pimgr[x-1+(y+2)*w+(z-1)*w*h], pimgr[x+(y+2)*w+(z-1)*w*h], pimgr[x+1+(y+2)*w+(z-1)*w*h], pimgr[x+2+(y+2)*w+(z-1)*w*h]};
                    type rx0[4] = {cubic(r00, dx), cubic(r01, dx), cubic(r02, dx), cubic(r03, dx)};

                    type r10[4] = {pimgr[x-1+(y-1)*w+z*w*h], pimgr[x+(y-1)*w+z*w*h], pimgr[x+1+(y-1)*w+z*w*h], pimgr[x+2+(y-1)*w+z*w*h]};
                    type r11[4] = {pimgr[x-1+(y)*w+z*w*h]  , pimgr[x+(y)*w+z*w*h]  , pimgr[x+1+(y)*w+z*w*h]  , pimgr[x+2+(y)*w+z*w*h]};
                    type r12[4] = {pimgr[x-1+(y+1)*w+z*w*h], pimgr[x+(y+1)*w+z*w*h], pimgr[x+1+(y+1)*w+z*w*h], pimgr[x+2+(y+1)*w+z*w*h]};
                    type r13[4] = {pimgr[x-1+(y+2)*w+z*w*h], pimgr[x+(y+2)*w+z*w*h], pimgr[x+1+(y+2)*w+z*w*h], pimgr[x+2+(y+2)*w+z*w*h]};
                    type rx1[4] = {cubic(r10, dx), cubic(r11, dx), cubic(r12, dx), cubic(r13, dx)};

                    type r20[4] = {pimgr[x-1+(y-1)*w+(z+1)*w*h], pimgr[x+(y-1)*w+(z+1)*w*h], pimgr[x+1+(y-1)*w+(z+1)*w*h], pimgr[x+2+(y-1)*w+(z+1)*w*h]};
                    type r21[4] = {pimgr[x-1+(y)*w+(z+1)*w*h]  , pimgr[x+(y)*w+(z+1)*w*h]  , pimgr[x+1+(y)*w+(z+1)*w*h]  , pimgr[x+2+(y)*w+(z+1)*w*h]};
                    type r22[4] = {pimgr[x-1+(y+1)*w+(z+1)*w*h], pimgr[x+(y+1)*w+(z+1)*w*h], pimgr[x+1+(y+1)*w+(z+1)*w*h], pimgr[x+2+(y+1)*w+(z+1)*w*h]};
                    type r23[4] = {pimgr[x-1+(y+2)*w+(z+1)*w*h], pimgr[x+(y+2)*w+(z+1)*w*h], pimgr[x+1+(y+2)*w+(z+1)*w*h], pimgr[x+2+(y+2)*w+(z+1)*w*h]};
                    type rx2[4] = {cubic(r20, dx), cubic(r21, dx), cubic(r22, dx), cubic(r23, dx)};

                    type r30[4] = {pimgr[x-1+(y-1)*w+(z+2)*w*h], pimgr[x+(y-1)*w+(z+2)*w*h], pimgr[x+1+(y-1)*w+(z+2)*w*h], pimgr[x+2+(y-1)*w+(z+2)*w*h]};
                    type r31[4] = {pimgr[x-1+(y)*w+(z+2)*w*h]  , pimgr[x+(y)*w+(z+2)*w*h]  , pimgr[x+1+(y)*w+(z+2)*w*h]  , pimgr[x+2+(y)*w+(z+2)*w*h]};
                    type r32[4] = {pimgr[x-1+(y+1)*w+(z+2)*w*h], pimgr[x+(y+1)*w+(z+2)*w*h], pimgr[x+1+(y+1)*w+(z+2)*w*h], pimgr[x+2+(y+1)*w+(z+2)*w*h]};
                    type r33[4] = {pimgr[x-1+(y+2)*w+(z+2)*w*h], pimgr[x+(y+2)*w+(z+2)*w*h], pimgr[x+1+(y+2)*w+(z+2)*w*h], pimgr[x+2+(y+2)*w+(z+2)*w*h]};
                    type rx3[4] = {cubic(r30, dx), cubic(r31, dx), cubic(r32, dx), cubic(r33, dx)};

                    type ry[4] = {cubic(rx0, dy), cubic(rx1, dy), cubic(rx2, dy), cubic(rx3, dy)};
                    pimgo[i + j*n0 + k*n0*n1] = cubic(ry, dz);
                }
                else if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l - 1) // linear interpolation otherwise
                {
                    type dx = xt - (type)x;
                    type dy = yt - (type)y;
                    type dxdy = dx*dy;
                    type rv = pimgr[x+y*w+z*w*h]*(1-dx-dy+dxdy) + pimgr[x+1+y*w+z*w*h]*(dx-dxdy) + pimgr[x+(y+1)*w+z*w*h]*(dy-dxdy) + pimgr[x+1+(y+1)*w+z*w*h]*dxdy;
                    type rw = pimgr[x+y*w+(z+1)*w*h]*(1-dx-dy+dxdy) + pimgr[x+1+y*w+(z+1)*w*h]*(dx-dxdy) + pimgr[x+(y+1)*w+(z+1)*w*h]*(dy-dxdy) + pimgr[x+1+(y+1)*w+(z+1)*w*h]*dxdy;
                    type dz = zt - (type)z;
                    type r = rv*(1-dz) + rw*dz;
                    pimgo[i + j*n0 + k*n0*n1] = r;
                }
                else if(x >= 0 && x < w && y >= 0 && y < h && z == l - 1) // border case
                {
                    type dx = xt - (type)x;
                    type dy = yt - (type)y;
                    type dxdy = dx*dy;
                    type rv = pimgr[x+y*w+z*w*h]*(1-dx-dy+dxdy) + pimgr[x+1+y*w+z*w*h]*(dx-dxdy) + pimgr[x+(y+1)*w+z*w*h]*(dy-dxdy) + pimgr[x+1+(y+1)*w+z*w*h]*dxdy;
                    pimgo[i + j*n0 + k*n0*n1] = rv;
                };
            };
        };
    };
};

#ifdef IMART_WITH_FFTW
template <typename type>
void vector_cpu<type>::fft(std::vector<pointer> & input, std::vector<pointer> & output, std::vector<int> size, bool forward)
{
    // fftw have limitations in size
    // fftw max 2d size 511*511
    int method;
    if (forward) method = FFTW_FORWARD;
    else method = FFTW_BACKWARD;

    // std::cout << "fft ";
    int dim = size.size();
    int N = 1;
    for(int i = 0; i < dim; i++) N *= size[i];
    
    // fftw_complex in[N], out[N];
    // Allocate memory with new for large vector size
    fftw_complex * in = new fftw_complex[N];
    fftw_complex * out = new fftw_complex[N];

    type * p1 = input[0]->data();
    type * p2 = input[1]->data();

    for (int i = 0; i < N; i++)
    {
        in[i][0] = p1[i];
        in[i][1] = p2[i];
    };

    #ifdef IMART_WITH_OPENMP
    if (N > 15*IMART_OPENMP_VECTOR_MIN_SIZE)
        fftw_plan_with_nthreads(omp_get_max_threads());
    #endif

    fftw_plan p_fft;
    // std::cout << "plan ";
    if (dim == 2) p_fft = fftw_plan_dft_2d(size[1], size[0], in, out, method, FFTW_ESTIMATE);
    if (dim == 3) p_fft = fftw_plan_dft_3d(size[2], size[1], size[0], in, out, method, FFTW_ESTIMATE);
    // std::cout << "plan ok ";
    fftw_execute(p_fft);
    fftw_destroy_plan(p_fft);

    type * pre = output[0]->data();
    type * pim = output[1]->data();
    if (forward)
    {
        for (int i = 0; i < N; i++)
        {
            pre[i] = out[i][0];
            pim[i] = out[i][1];
        };
    }
    else // if ifft divide by N
    {
        for (int i = 0; i < N; i++)
        {
            pre[i] = out[i][0];
            pim[i] = out[i][1];
        };
    }

    delete[] in;
    delete[] out;
};
#endif

template <typename type>
void vector_cpu<type>::gradientx(typename vector_cpu<type>::pointer imgr,
                                 typename vector_cpu<type>::pointer imgo,
                                 std::vector<int> ref_size)
{
    type a = 0.5; // finite difference scalar
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    if (ref_size.size() == 2)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1];
        #pragma omp parallel for collapse(2) if (n0*n1 > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int j = 0; j < n1; j++)
        {
            for(int i = 0; i < n0; i++)
            {
                if(i == 0)          pimgo[i + j*n0] = pimgr[i+1 + j*n0] - pimgr[i + j*n0];
                else if (i == n0-1) pimgo[i + j*n0] = pimgr[i + j*n0] - pimgr[i-1 + j*n0];
                else      pimgo[i + j*n0] = a*pimgr[i+1 + j*n0] - a*pimgr[i-1 + j*n0];
            };
        };
    }
    else if (ref_size.size() == 3)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[2];
        #pragma omp parallel for collapse(3) if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int k = 0; k < n2; k++)
        {
            for(int j = 0; j < n1; j++)
            {
                for(int i = 0; i < n0; i++)
                {
                    if(i == 0)          pimgo[i + j*n0 + k*n0*n1] = pimgr[i+1 + j*n0 + k*n0*n1] - pimgr[i + j*n0 + k*n0*n1];
                    else if (i == n0-1) pimgo[i + j*n0 + k*n0*n1] = pimgr[i + j*n0 + k*n0*n1] - pimgr[i-1 + j*n0 + k*n0*n1];
                    else      pimgo[i + j*n0 + k*n0*n1] = a*pimgr[i+1 + j*n0 + k*n0*n1] - a*pimgr[i-1 + j*n0 + k*n0*n1];
                };
            };
        };
    };
};

template <typename type>
void vector_cpu<type>::gradienty(typename vector_cpu<type>::pointer imgr,
                                 typename vector_cpu<type>::pointer imgo,
                                 std::vector<int> ref_size)
{
    type a = 0.5; // finite difference scalar
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    if (ref_size.size() == 2)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1];
        #pragma omp parallel for collapse(2) if (n0*n1 > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int j = 0; j < n1; j++)
        {
            for(int i = 0; i < n0; i++)
            {
                if(j == 0)          pimgo[i + j*n0] = pimgr[i + (j+1)*n0] - pimgr[i + j*n0];
                else if (j == n1-1) pimgo[i + j*n0] = pimgr[i + j*n0] - pimgr[i + (j-1)*n0];
                else    pimgo[i + j*n0] = a*pimgr[i + (j+1)*n0] - a*pimgr[i + (j-1)*n0];
            };
        };
    }
    else if (ref_size.size() == 3)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[2];
        #pragma omp parallel for collapse(3) if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int k = 0; k < n2; k++)
        {
            for(int j = 0; j < n1; j++)
            {
                for(int i = 0; i < n0; i++)
                {
                    if(j == 0)          pimgo[i + j*n0 + k*n0*n1] = pimgr[i + (j+1)*n0 + k*n0*n1] - pimgr[i + j*n0 + k*n0*n1];
                    else if (j == n1-1) pimgo[i + j*n0 + k*n0*n1] = pimgr[i + j*n0 + k*n0*n1] - pimgr[i + (j-1)*n0 + k*n0*n1];
                    else    pimgo[i + j*n0 + k*n0*n1] = a*pimgr[i + (j+1)*n0 + k*n0*n1] - a*pimgr[i + (j-1)*n0 + k*n0*n1];
                };
            };
        };
    };
};

template <typename type>
void vector_cpu<type>::gradientz(typename vector_cpu<type>::pointer imgr,
                                 typename vector_cpu<type>::pointer imgo,
                                 std::vector<int> ref_size)
{
    type a = 0.5; // finite difference scalar
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    if (ref_size.size() == 3)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[2];
        #pragma omp parallel for collapse(3) if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int k = 0; k < n2; k++)
        {
            for(int j = 0; j < n1; j++)
            {
                for(int i = 0; i < n0; i++)
                {
                    if(k == 0)          pimgo[i + j*n0 + k*n0*n1] = pimgr[i + j*n0 + (k+1)*n0*n1] - pimgr[i + j*n0 + k*n0*n1];
                    else if (k == n2-1) pimgo[i + j*n0 + k*n0*n1] = pimgr[i + j*n0 + k*n0*n1] - pimgr[i + j*n0 + (k-1)*n0*n1];
                    else    pimgo[i + j*n0 + k*n0*n1] = a*pimgr[i + j*n0 + (k+1)*n0*n1] - a*pimgr[i + j*n0 + (k-1)*n0*n1];
                };
            };
        };
    };
};

template <typename type>
void vector_cpu<type>::convolution( typename vector_cpu<type>::pointer imgr,
                                    typename vector_cpu<type>::pointer kernel,
                                    typename vector_cpu<type>::pointer imgo,
                                    std::vector<int> ref_size, std::vector<int> kernel_size)
{
    type * pimgr = imgr->data();
    type * pkrnl = kernel->data();
    type * pimgo = imgo->data();

    if (ref_size.size() == 2)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1];
        int kw0 = kernel_size[0]; int kw1 = kernel_size[1];
        int off0 = kw0>>1; int off1 = kw1>>1;

        #pragma omp parallel for collapse(2) if (n0*n1 > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int j = off1; j < n1-off1; j++)
        {
            for(int i = off0; i < n0-off0; i++)
            {
                type sum = 0;
                for (int p = 0; p < kw1; p++)
                {
                    for (int q = 0; q < kw0; q++)
                    {
                        sum += pimgr[i+q-off0 + (j+p-off1)*n0] * pkrnl[p*kw0 + q];
                    };
                };
                pimgo[i + j*n0] = sum;
            };
        };
    }
    else if (ref_size.size() == 3)
    {
        // printf("n0: %4d, n1: %4d, n2: %4d\n", n0, n1, n2);
        int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[2];
        int kw0 = kernel_size[0]; int kw1 = kernel_size[1]; int kw2 = kernel_size[2];
        int off0 = kw0>>1; int off1 = kw1>>1; int off2 = kw2>>1;
        #pragma omp parallel for collapse(3) if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int k = off2; k < n2-off2; k++)
        {
            for(int j = off1; j < n1-off1; j++)
            {
                for(int i = off0; i < n0-off0; i++)
                {
                    type sum = 0;
                    for (int r = 0; r < kw2; r++)
                    {
                        for (int p = 0; p < kw1; p++)
                        {
                            for (int q = 0; q < kw0; q++)
                            {
                                sum += pimgr[i+q-off0 + (j+p-off1)*n0 + (k+r-off2)*n0*n1] * pkrnl[r*kw0*kw1 + p*kw0 + q];
                            };
                        };
                    };
                    pimgo[i + j*n0 + k*n0*n1] = sum;
                };
            };
        };
    };
};

/*
template <typename type>
void vector_cpu<type>::convolution( typename vector_cpu<type>::pointer imgr,
                                    typename vector_cpu<type>::pointer kernel,
                                    typename vector_cpu<type>::pointer imgo,
                                    std::vector<int> ref_size, std::vector<int> kernel_size)
{
    type * pimgr = imgr->data();
    type * pkrnl = kernel->data();
    type * pimgo = imgo->data();

    if (ref_size.size() == 2)
    {
        int dataSizeX = ref_size[0]; int dataSizeY = ref_size[1];
        int kernelSizeX = kernel_size[0]; int kernelSizeY = kernel_size[1];
        
        type *inPtr, *inPtr2, *outPtr, *kPtr;

        int kCenterX, kCenterY;
        int rowMin, rowMax;                             // to check boundary of input array
        int colMin, colMax;                             //
        
        // find center position of kernel (half of kernel size)
        kCenterX = kernelSizeX >> 1;
        kCenterY = kernelSizeY >> 1;

        // init working  pointers
        inPtr = inPtr2 = &pimgr[dataSizeX * kCenterY + kCenterX];  // note that  it is shifted (kCenterX, kCenterY),
        outPtr = pimgo;
        kPtr = pkrnl;

        // start convolution
        #pragma omp parallel for collapse(2) if (dataSizeX*dataSizeY > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int i = 0; i < dataSizeY; ++i)                   // number of rows
        {
            for(int j = 0; j < dataSizeX; ++j)              // number of columns
            {
                // compute the range of convolution, the current row of kernel should be between these
                rowMax = i + kCenterY;
                rowMin = i - dataSizeY + kCenterY;
                
                // compute the range of convolution, the current column of kernel should be between these
                colMax = j + kCenterX;
                colMin = j - dataSizeX + kCenterX;

                *outPtr = 0;                            // set to 0 before accumulate

                // flip the kernel and traverse all the kernel values
                // multiply each kernel value with underlying input data
                for(int m = 0; m < kernelSizeY; ++m)        // kernel rows
                {
                    // check if the index is out of bound of input array
                    if(m <= rowMax && m > rowMin)
                    {
                        for(int n = 0; n < kernelSizeX; ++n)
                        {
                            // check the boundary of array
                            if(n <= colMax && n > colMin)
                                *outPtr += *(inPtr - n) * *kPtr;
                            ++kPtr;                     // next kernel
                        }
                    }
                    else
                        kPtr += kernelSizeX;            // out of bound, move to next row of kernel

                    inPtr -= dataSizeX;                 // move input data 1 raw up
                }

                kPtr = pkrnl;                          // reset kernel to (0,0)
                inPtr = ++inPtr2;                       // next input
                ++outPtr;                               // next output
            }
        };
    }
    else if (ref_size.size() == 3)
    {
        // printf("n0: %4d, n1: %4d, n2: %4d\n", n0, n1, n2);
        int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[2];
        int kw0 = kernel_size[0]; int kw1 = kernel_size[1]; int kw2 = kernel_size[2];
        int off0 = kw0>>1; int off1 = kw1>>1; int off2 = kw2>>1;
        #pragma omp parallel for collapse(3) if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int k = off2; k < n2-off2; k++)
        {
            for(int j = off1; j < n1-off1; j++)
            {
                for(int i = off0; i < n0-off0; i++)
                {
                    type sum = 0;
                    for (int r = 0; r < kw2; r++)
                    {
                        for (int p = 0; p < kw1; p++)
                        {
                            for (int q = 0; q < kw0; q++)
                            {
                                sum += pimgr[i+q-off0 + (j+p-off1)*n0 + (k+r-off2)*n0*n1] * pkrnl[r*kw0*kw1 + p*kw0 + q];
                            };
                        };
                    };
                    pimgo[i + j*n0 + k*n0*n1] = sum;
                };
            };
        };
    };
};*/

// template <typename type>
// void vector_cpu<type>::convolution( typename vector_cpu<type>::pointer imgr,
//                                     typename vector_cpu<type>::pointer kernel,
//                                     typename vector_cpu<type>::pointer imgo,
//                                     std::vector<int> ref_size, int kwidth)
// {
//     type * pimgr = imgr->data();
//     type * pkrnl = kernel->data();
//     type * pimgo = imgo->data();

//     int off = std::floor(kwidth/2.0);

//     if (ref_size.size() == 2)
//     {
//         int n0 = ref_size[0]; int n1 = ref_size[1];
//         #pragma omp parallel for collapse(2) if (n0*n1 > IMART_OPENMP_VECTOR_MIN_SIZE)
//         for(int j = off; j < n1-off; j++)
//         {
//             for(int i = off; i < n0-off; i++)
//             {
//                 type sum = 0;
//                 for (int p = 0; p < kwidth; p++)
//                 {
//                     for (int q = 0; q < kwidth; q++)
//                     {
//                         sum += pimgr[i+q-off + (j+p-off)*n0] * pkrnl[p*kwidth + q];
//                     };
//                 };
//                 pimgo[i + j*n0] = sum;
//             };
//         };
//     }
//     else if (ref_size.size() == 3)
//     {
//         // printf("n0: %4d, n1: %4d, n2: %4d\n", n0, n1, n2);
//         int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[2];
//         #pragma omp parallel for collapse(3) if (n0*n1*n2 > IMART_OPENMP_VECTOR_MIN_SIZE)
//         for(int k = off; k < n2-off; k++)
//         {
//             for(int j = off; j < n1-off; j++)
//             {
//                 for(int i = off; i < n0-off; i++)
//                 {
//                     type sum = 0;
//                     for (int r = 0; r < kwidth; r++)
//                     {
//                         for (int p = 0; p < kwidth; p++)
//                         {
//                             for (int q = 0; q < kwidth; q++)
//                             {
//                                 sum += pimgr[i+q-off + (j+p-off)*n0 + (k+r-off)*n0*n1] * pkrnl[r*kwidth*kwidth + p*kwidth + q];
//                             };
//                         };
//                     };
//                     pimgo[i + j*n0 + k*n0*n1] = sum;
//                 };
//             };
//         };
//     };
// };

}; //end namespace

#endif