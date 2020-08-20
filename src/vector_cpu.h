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

// fftw
#include <fftw3.h>      // fft library

// local libs
#include "inherit.h"
#include "object.h"

namespace imart
{

// Class object
template <typename type>
class vector_cpu: public inherit<vector_cpu<type>, object>, std::vector<type>
{
public:
    //Type definitions
    using self    = vector_cpu;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using object::class_name;
    using object::get_name;
    using std::vector<type>::size;
    using std::vector<type>::begin;
    using std::vector<type>::end;
    using std::vector<type>::operator[];
    using std::vector<type>::data;

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
    vector_cpu(): std::vector<type>() { class_name = "vector_cpu"; };          // constructor empty
    vector_cpu(int s): std::vector<type>(s) { class_name = "vector_cpu"; };    // constructor
    vector_cpu(int s, type value): std::vector<type>(s, value) { class_name = "vector_cpu"; };
    vector_cpu(std::initializer_list<type> list): std::vector<type>(list) { class_name = "vector_cpu"; };
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

    static void fft(std::vector<pointer> & input, std::vector<pointer> & output, std::vector<int> size, bool forward);

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
                             std::vector<int> ref_size, int kwidth);

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
    const type * p2 = input.data();

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
    return static_cast<std::vector<type>>(*this);
};

// ===========================================
// Memory Functions
// ===========================================
template <typename type>
void vector_cpu<type>::read_ram(type * p, int s, int offset)
{
    type * d = this->data() + offset;
    for(int k=0; k<s; k++)
    {
        *(d+k) = *(p+k);
    };
};

template <typename type>
void vector_cpu<type>::write_ram(type * p, int s, int offset)
{
    type * d = this->data() + offset;
    for(int k=0; k<s; k++)
    {
        *(p+k) = *(d+k);
    };
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
    assert(this->size() == input.size());
};

// ===========================================
// Initialization Functions
// ===========================================
template <typename type>
void vector_cpu<type>::zeros()
{
    int size = this->size();
    type * p = this->data();
    for(int k=0; k<size; k++)
    {
        p[k] = (type)0; // casting to pixel_type
    };
};

template <typename type>
void vector_cpu<type>::ones()
{
    int size = this->size();
    type * p = this->data();
    for(int k=0; k<size; k++)
    {
        p[k] = (type)1.0; // casting to pixel_type
    };
};

// Set all pixel to a fixed value
template <typename type>
void vector_cpu<type>::assign(type value)
{
    int size = this->size();
    type * p = this->data();
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
    for(int k=0; k<size; k++)
    {
        p[k] = (type)uniform(gen); // casting to pixel_type
    };
};

// ===========================================
// Overloading operators
// ===========================================
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
    const type * p2 = input.data();
    type * p3 = output->data();

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
    const type * p2 = input.data();
    type * p3 = output->data();

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
    const type * p2 = input.data();
    type * p3 = output->data();

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
    const type * p2 = input.data();
    type * p3 = output->data();

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
    const type * p2 = input.data();
    type * p3 = output->data();

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
    const type * p1 = input.data();
    type * p2 = output->data();

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
    const type * p1 = input.data();
    type * p2 = output->data();

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
    const type * p1 = input.data();
    type * p2 = output->data();

    for(int k=0; k<size; k++)
    {
        p2[k] = pow(scalar, p1[k]);
    };
    return output;
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
    const type * p2 = input.data();
    int size = this->size();

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

    // #pragma omp parallel for
    int sz = xo->size();
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

    // #pragma omp parallel for
    int sz = xo->size();
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

    // #pragma omp parallel for
    int sz = xo->size();
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

    // #pragma omp parallel for
    int sz = xo->size();
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
    // #pragma omp parallel for
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
    // #pragma omp parallel for
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

    // #pragma omp parallel for
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

    // #pragma omp parallel for
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

    for(int j = 0; j < n1; j++)
    {
        for(int i = 0; i < n0; i++)
        {
            type xt = pxo[i + j*n0];
            type yt = pyo[i + j*n0];
            int x = floor(xt);
            int y = floor(yt);

            if(x >= 0 && x < w && y >= 0 && y < h - 1)
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
                if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l - 1)
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
    
    fftw_complex in[N], out[N];
    type * p1 = input[0]->data();
    type * p2 = input[1]->data();

    for (int i = 0; i < N; i++)
    {
        in[i][0] = p1[i];
        in[i][1] = p2[i];
    };

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
            pre[i] = out[i][0]/N;
            pim[i] = out[i][1]/N;
        };
    }
};

template <typename type>
void vector_cpu<type>::gradientx(typename vector_cpu<type>::pointer imgr,
                                 typename vector_cpu<type>::pointer imgo,
                                 std::vector<int> ref_size)
{
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    if (ref_size.size() == 2)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1];
        for(int j = 0; j < n1; j++)
        {
            for(int i = 0; i < n0; i++)
            {
                if(i == 0)          pimgo[i + j*n0] = pimgr[i+1 + j*n0] - pimgr[i + j*n0];
                else if (i == n0-1) pimgo[i + j*n0] = pimgr[i + j*n0] - pimgr[i-1 + j*n0];
                else      pimgo[i + j*n0] = 0.5*pimgr[i+1 + j*n0] - 0.5*pimgr[i-1 + j*n0];
            };
        };
    }
    else if (ref_size.size() == 3)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[1];
        for(int k = 0; k < n2; k++)
        {
            for(int j = 0; j < n1; j++)
            {
                for(int i = 0; i < n0; i++)
                {
                    if(i == 0)          pimgo[i + j*n0 + k*n0*n1] = pimgr[i+1 + j*n0 + k*n0*n1] - pimgr[i + j*n0 + k*n0*n1];
                    else if (i == n0-1) pimgo[i + j*n0 + k*n0*n1] = pimgr[i + j*n0 + k*n0*n1] - pimgr[i-1 + j*n0 + k*n0*n1];
                    else      pimgo[i + j*n0 + k*n0*n1] = 0.5*pimgr[i+1 + j*n0 + k*n0*n1] - 0.5*pimgr[i-1 + j*n0 + k*n0*n1];
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
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    if (ref_size.size() == 2)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1];
        for(int j = 0; j < n1; j++)
        {
            for(int i = 0; i < n0; i++)
            {
                if(j == 0)          pimgo[i + j*n0] = pimgr[i + (j+1)*n0] - pimgr[i + j*n0];
                    else if (j == n1-1) pimgo[i + j*n0] = pimgr[i + j*n0] - pimgr[i + (j-1)*n0];
                    else    pimgo[i + j*n0] = 0.5*pimgr[i + (j+1)*n0] - 0.5*pimgr[i + (j-1)*n0];
            };
        };
    }
    else if (ref_size.size() == 3)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[1];
        for(int k = 0; k < n2; k++)
        {
            for(int j = 0; j < n1; j++)
            {
                for(int i = 0; i < n0; i++)
                {
                    if(j == 0)          pimgo[i + j*n0 + k*n0*n1] = pimgr[i + (j+1)*n0 + k*n0*n1] - pimgr[i + j*n0 + k*n0*n1];
                    else if (j == n1-1) pimgo[i + j*n0 + k*n0*n1] = pimgr[i + j*n0 + k*n0*n1] - pimgr[i + (j-1)*n0 + k*n0*n1];
                    else    pimgo[i + j*n0 + k*n0*n1] = 0.5*pimgr[i + (j+1)*n0 + k*n0*n1] - 0.5*pimgr[i + (j-1)*n0 + k*n0*n1];
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
    type * pimgr = imgr->data();
    type * pimgo = imgo->data();

    if (ref_size.size() == 3)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[1];
        for(int k = 0; k < n2; k++)
        {
            for(int j = 0; j < n1; j++)
            {
                for(int i = 0; i < n0; i++)
                {
                    if(k == 0)          pimgo[i + j*n0 + k*n0*n1] = pimgr[i + j*n0 + (k+1)*n0*n1] - pimgr[i + j*n0 + k*n0*n1];
                    else if (k == n2-1) pimgo[i + j*n0 + k*n0*n1] = pimgr[i + j*n0 + k*n0*n1] - pimgr[i + j*n0 + (k-1)*n0*n1];
                    else    pimgo[i + j*n0 + k*n0*n1] = 0.5*pimgr[i + j*n0 + (k+1)*n0*n1] - 0.5*pimgr[i + j*n0 + (k-1)*n0*n1];
                };
            };
        };
    };
};

template <typename type>
void vector_cpu<type>::convolution( typename vector_cpu<type>::pointer imgr,
                                    typename vector_cpu<type>::pointer kernel,
                                    typename vector_cpu<type>::pointer imgo,
                                    std::vector<int> ref_size, int kwidth)
{
    type * pimgr = imgr->data();
    type * pkrnl = kernel->data();
    type * pimgo = imgo->data();

    int off = std::floor(kwidth/2.0);

    if (ref_size.size() == 2)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1];
        
        for(int j = off; j < n1-off; j++)
        {
            for(int i = off; i < n0-off; i++)
            {
                type sum = 0;
                for (int p = 0; p < kwidth; p++)
                {
                    for (int q = 0; q < kwidth; q++)
                    {
                        sum += pimgr[i+p-off + (j+q-off)*n0] * pkrnl[p*kwidth + q];
                    };
                };
                pimgo[i + j*n0] = sum;
            };
        };
    }
    else if (ref_size.size() == 3)
    {
        int n0 = ref_size[0]; int n1 = ref_size[1]; int n2 = ref_size[1];
        
        for(int k = off; k < n2-off; k++)
        {
            for(int j = off; j < n1-off; j++)
            {
                for(int i = off; i < n0-off; i++)
                {
                    type sum = 0;
                    for (int r = 0; r < kwidth; r++)
                    {
                        for (int p = 0; p < kwidth; p++)
                        {
                            for (int q = 0; q < kwidth; q++)
                            {
                                sum += pimgo[i+p-off + (j+q-off)*n0 + (k-off)*n0*n1] * pkrnl[r*kwidth*kwidth + p*kwidth + q];
                            };
                        };
                    };
                    pimgo[i + j*n0 + k*n0*n1] = sum;
                };
            };
        };
    };
};

}; //end namespace

#endif