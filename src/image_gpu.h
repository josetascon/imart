/*
* @Author: jose
* @Date:   2020-03-03 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-03 00:00:00
*/

#ifndef __IMAGE_GPU_H__
#define __IMAGE_GPU_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <random>       // std::random

// gpu libs
#include <CL/cl.hpp>

// ViennaCL headers
#include <viennacl/scalar.hpp>
#include <viennacl/matrix.hpp>

// images itk
#include <itkImage.h>
#include <itkImageFileReader.h>

// local libs
#include "object.h"
#include "image_base.h"

namespace imart
{

// Class image_gpu
template <typename pixel_type>
class image_gpu: public image_base<pixel_type>
{
public:
    //Type definitions
    using iterator       = typename viennacl::const_vector_iterator<pixel_type, 1>;
    using ptr_vector     = std::shared_ptr<std::vector<pixel_type>>;
    using ptr_vcl_vector = std::shared_ptr< viennacl::vector<pixel_type> >;
    // using ptr_pixels4 = std::unique_ptr<std::array<pixel_type,4>>;
    // using ptr_pixels8 = std::unique_ptr<std::array<pixel_type,8>>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    using image_base<pixel_type>::width;            // consider change image_base as virtual only, image_base cpu methods move to cpu
    using image_base<pixel_type>::height;
    using image_base<pixel_type>::length;
    using image_base<pixel_type>::num_elements;
    using image_base<pixel_type>::channels;         // for now will only support one channel
    
    // Image data
    ptr_vcl_vector data;
    // MISSIN CODE. TAKE IT FROM image_base

    void allocate(int elements);
    void init(int w, int h);
    

public:
    //
    image_gpu();
    //! Constructor 
    image_gpu(int w, int h);        // Only 2d
    //! Constructor with existing data.
    image_gpu(ptr_vector buffer, int w, int h);        // Only 2d

    image_gpu(const image_gpu<pixel_type> & input);

    //! Full copy image
    void copy(const image_gpu & input);
    //! Duplicate data array for images
    void duplicate(const image_gpu & input);
    //! Imitate other image with a property copy
    void imitate(const image_gpu & input);

    using image_base<pixel_type>::assert_size;

    std::string info(std::string msg);
    std::string info_data(std::string msg);

    int get_ptr_count() const;
    ptr_vcl_vector get_data() const;

    iterator begin() const;

    iterator end() const;

    // void assert_size(const image_gpu<pixel_type> & input);

    // ===========================================
    // Initialization Functions
    // ===========================================
    void zeros();
    void ones();
    void identity(); // TODO
    void random(float min=0.0, float max=1.0);
    void fill(pixel_type value);

    // ===========================================
    // Overloading Operators
    // ===========================================
    // Operators
    image_gpu<pixel_type> & operator = (const image_gpu<pixel_type> & input);
    
    // Image to Image
    image_gpu<pixel_type> operator + (const image_gpu<pixel_type> & input);
    image_gpu<pixel_type> operator - (const image_gpu<pixel_type> & input);
    image_gpu<pixel_type> operator * (const image_gpu<pixel_type> & input);
    image_gpu<pixel_type> operator / (const image_gpu<pixel_type> & input);
    image_gpu<pixel_type> operator ^ (const image_gpu<pixel_type> & input);

    // Scalar
    image_gpu<pixel_type> operator + (pixel_type scalar);
    image_gpu<pixel_type> operator - (pixel_type scalar);
    image_gpu<pixel_type> operator * (pixel_type scalar);
    image_gpu<pixel_type> operator / (pixel_type scalar);
    image_gpu<pixel_type> operator ^ (pixel_type scalar);

    // Friend classes to support double side
    template<typename pixel_t>
    friend image_gpu<pixel_t> operator + (pixel_t scalar, image_gpu<pixel_t> & input);
    template<typename pixel_t>
    friend image_gpu<pixel_t> operator - (pixel_t scalar, const image_gpu<pixel_t> & input);
    template<typename pixel_t>
    friend image_gpu<pixel_t> operator * (pixel_t scalar, image_gpu<pixel_t> & input);
    template<typename pixel_t>
    friend image_gpu<pixel_t> operator / (pixel_t scalar, const image_gpu<pixel_t> & input);

    // Reduction operations
    pixel_type min();
    pixel_type max();
    pixel_type sum();
    pixel_type prod();


};

template <typename pixel_type>
image_gpu<pixel_type>::image_gpu()
{
    init( 0, 0 );
    allocate(0);
};

template <typename pixel_type>
image_gpu<pixel_type>::image_gpu(int w, int h)
{
    init( w, h );
    allocate(w*h);
};

template <typename pixel_type>
image_gpu<pixel_type>::image_gpu(ptr_vector buffer, int w, int h)
{
    init( w, h );
    allocate( w*h );
    // std::cout << "allocate" << std::endl;
    viennacl::copy(buffer->begin(), buffer->end(), data->begin());
    // std::cout << "copy" << std::endl;
};

// Constructor
template <typename pixel_type>
image_gpu<pixel_type>::image_gpu(const image_gpu<pixel_type> & input)
{
    copy(input);
};

template <typename pixel_type>
void image_gpu<pixel_type>::allocate(int elements)
{
    data = std::make_shared<viennacl::vector<pixel_type>>(elements);
    // std::cout << data;
};

template <typename pixel_type>
void image_gpu<pixel_type>::copy(const image_gpu<pixel_type> & input)
{
    this->copy_properties(input);
    allocate(input.get_total_elements());
    viennacl::copy(input.begin(), input.end(), data->begin());
};

// Point to the same data
template <typename pixel_type>
void image_gpu<pixel_type>::duplicate(const image_gpu<pixel_type> & input)
{
    this->copy_properties(input);
    data.reset();
    data = input.get_data();
};

template <typename pixel_type>
void image_gpu<pixel_type>::imitate(const image_gpu<pixel_type> & input)
{
    this->copy_properties(input);
    data.reset();
    allocate(input.get_total_elements());
};

template <typename pixel_type>
void image_gpu<pixel_type>::init(int w, int h)
{   
    this->class_name = "image_gpu";
    channels = 1;
    width = w;
    height = h;
    length = 1;
    num_elements = width*height;

    object<pixel_type>::init(2);

    this->size = std::vector<int>{width, height};
};


template <typename pixel_type>
std::string image_gpu<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Image Information";
    if (msg != "") { title = msg; };
    // Summary of the image information
    ss << object<pixel_type>::info(title);
    // ss << "\n===== " << title << " =====\n";
    // ss << "Pixel type: \t\t" << typeid((*data)[0]).name() << std::endl;
    // ss << "Dimensions: \t\t" << dim << std::endl;
    ss << "Pixel channels: \t" << channels << std::endl; 
    
    ss << "Size: \t\t\t[ ";
    for(int i = 0; i < this->size.size(); i++) { ss << this->size[i] << " "; };
    ss << "]" << std::endl;
    
    ss << "Length (mm): \t\t[ ";
    for(int i = 0; i < this->spacing.size(); i++) { ss << this->spacing[i] << " "; };
    ss << "]" << std::endl;

    ss << "Origin (mm): \t\t[ ";
    for(int i = 0; i < this->origin.size(); i++) { ss << this->origin[i] << " "; };
    ss << "]" << std::endl;

    ss << "Total elements: \t" << num_elements << std::endl; //Get the total number of pixels

    ss << "Data pointer: \t\t" << get_data() << std::endl; //Get the total number of pixels
    // ss << std::endl;

    return ss.str();
};

template <typename pixel_type>
std::string image_gpu<pixel_type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; };
    ss << *data;
    ss << std::endl;
    return ss.str();
};

template <typename pixel_type>
typename image_gpu<pixel_type>::ptr_vcl_vector image_gpu<pixel_type>::get_data() const
{
    return data;
};

template <typename pixel_type>
int image_gpu<pixel_type>::get_ptr_count() const
{
    return data.use_count();
};

template <typename pixel_type>
typename image_gpu<pixel_type>::iterator image_gpu<pixel_type>::begin() const
{
    return data.get()->begin();
};

template <typename pixel_type>
typename image_gpu<pixel_type>::iterator image_gpu<pixel_type>::end() const
{
    return data.get()->end();
};

// ===========================================
// Initialization Functions
// ===========================================
template <typename pixel_type>
void image_gpu<pixel_type>::zeros()
{
    *data = viennacl::zero_vector<pixel_type>(num_elements);
};

template <typename pixel_type>
void image_gpu<pixel_type>::ones()
{
    *data = viennacl::scalar_vector<pixel_type>(num_elements, 1);
};

template <typename pixel_type>
void image_gpu<pixel_type>::random(float min, float max)
{
    // viennacl::tools::uniform_random_numbers<pixel_type> randomNumber;
    // *data = viennacl::random_vector<pixel_type>(num_elements, randomNumber);
};

// Set all pixel to a fixed value
template <typename pixel_type>
void image_gpu<pixel_type>::fill(pixel_type value)
{
    *data = viennacl::scalar_vector<pixel_type>(num_elements, value);
};

// Equal
template <typename pixel_type>
image_gpu<pixel_type> & image_gpu<pixel_type>::operator = (const image_gpu<pixel_type> & input)
{
    // delete &data;
    duplicate(input);
    return *this;
};

// Image to Image
template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator + (const image_gpu<pixel_type> & input)
{
    assert_size(input);

    image_gpu<pixel_type> result;
    result.imitate(input);
    *(result.get_data()) = *(this->get_data()) + *(input.get_data());
    
    return result;
};

template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator - (const image_gpu<pixel_type> & input)
{
    assert_size(input);

    image_gpu<pixel_type> result;
    result.imitate(input);
    *(result.get_data()) = *(this->get_data()) - *(input.get_data());
    
    return result;
};

template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator * (const image_gpu<pixel_type> & input)
{
    assert_size(input);

    image_gpu<pixel_type> result;
    result.imitate(input);
    *(result.get_data()) = *(this->get_data()) * *(input.get_data());
    
    return result;
};

template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator / (const image_gpu<pixel_type> & input)
{
    assert_size(input);

    image_gpu<pixel_type> result;
    result.imitate(input);
    *(result.get_data()) = viennacl::linalg::element_div( *(this->get_data()) , *(input.get_data()) );
    
    return result;
};

template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator ^ (const image_gpu<pixel_type> & input)
{
    assert_size(input);

    image_gpu<pixel_type> result;
    result.imitate(input);

    *(result.get_data()) = viennacl::linalg::element_pow(*(this->get_data()), *(input.get_data()));
    
    return result;
};

template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator + (pixel_type scalar)
{
    image_gpu<pixel_type> result;
    result.imitate(*this);
    result.fill(scalar);
    // viennacl::scalar<pixel_type> ss = scalar;
    // *(result.get_data()) = *(this->get_data()) + ss;
    // *(result.get_data()) += ss;

    result = result + *this;

    return result;
};

template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator - (pixel_type scalar)
{
    image_gpu<pixel_type> result;
    result.imitate(*this);
    result.fill(scalar);
    // *(result.get_data()) = *(this->get_data()) - scalar;

    result = *this - result;

    return result;
};

template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator * (pixel_type scalar)
{
    image_gpu<pixel_type> result;
    result.imitate(*this);
    viennacl::scalar<pixel_type> ss = scalar;
    *(result.get_data()) = *(this->get_data()) * ss;

    return result;
};

template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator / (pixel_type scalar)
{
    image_gpu<pixel_type> result;
    result.imitate(*this);
    viennacl::scalar<pixel_type> ss = 1/scalar;
    *(result.get_data()) = *(this->get_data()) * ss;

    return result;
};

template <typename pixel_type>
image_gpu<pixel_type> image_gpu<pixel_type>::operator ^ (pixel_type scalar)
{
    image_gpu<pixel_type> result;
    result.imitate(*this);
    *(result.get_data()) = viennacl::linalg::element_pow(*(this->get_data()), scalar);

    return result;
};

// Scalar left hand side defined in header due to use of friend functions
template <typename pixel_type>
image_gpu<pixel_type> operator + (pixel_type scalar, image_gpu<pixel_type> & input)
{
    return input + scalar;
};

template <typename pixel_type>
image_gpu<pixel_type> operator - (pixel_type scalar, const image_gpu<pixel_type> & input)
{
    image_gpu<pixel_type> result;
    result.imitate(input);
    result.fill(scalar);
    return result - input;
};

template <typename pixel_type>
image_gpu<pixel_type> operator * (pixel_type scalar, image_gpu<pixel_type> & input)
{
    return input * scalar;
};

template <typename pixel_type>
image_gpu<pixel_type> operator / (pixel_type scalar, const image_gpu<pixel_type> & input)
{
    // return input * scalar;
    image_gpu<pixel_type> result;
    result.imitate(input);
    result.fill(scalar);
    return result / input;
};

// Reductions operations
// CONTINUE HERE
template <typename pixel_type>
pixel_type image_gpu<pixel_type>::min()
{
    pixel_type x = 0;
    return x;
};

template <typename pixel_type>
pixel_type image_gpu<pixel_type>::max()
{
    pixel_type x = 0;
    
    return x;
};

template <typename pixel_type>
pixel_type image_gpu<pixel_type>::sum()
{
    pixel_type x = 0;
    
    return x;
};

template <typename pixel_type>
pixel_type image_gpu<pixel_type>::prod()
{
    pixel_type x = 1;
    
    return x;
};

}; //end namespace

#endif