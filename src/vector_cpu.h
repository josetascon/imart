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

    static std::vector<typename vector_cpu<type>::pointer> grid_2d(int w, int h, std::vector<double> & sod);
    static std::vector<typename vector_cpu<type>::pointer> grid_3d(int w, int h, int l, std::vector<double> & sod);
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
void vector_cpu<type>::read_ram(type * p, int s)
{
    type * d = this->data();
    for(int k=0; k<s; k++)
    {
        *(d+k) = *(p+k);
    };
};

template <typename type>
void vector_cpu<type>::write_ram(type * p, int s)
{
    type * d = this->data();
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
std::vector<typename vector_cpu<type>::pointer> vector_cpu<type>::grid_2d(int w, int h, std::vector<double> & sod)
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
std::vector<typename vector_cpu<type>::pointer> vector_cpu<type>::grid_3d(int w, int h, int l, std::vector<double> & sod)
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


}; //end namespace

#endif