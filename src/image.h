/*
* @Author: jose
* @Date:   2019-11-05 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-14 17:39:46
*/

#ifndef __IMAGE_H__
#define __IMAGE_H__

// std libs
#include <iostream>     // std::cout
#include <sstream>      // stringstream
#include <memory>       // smart pointers
#include <vector>       // std::vector

// itk headers
#include <itkPoint.h>
#include <itkImage.h>
#include <itkImageBase.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// extra matrix eigen
// #include <eigen3/Eigen/Core>

// local libs
#include "space_object.h"
#include "data_object.h"
#include "vector_cpu.h"
#ifdef IMART_WITH_OPENCL
#include "vector_opencl.h"
#endif
#ifdef IMART_WITH_CUDA
#include "vector_cuda.h"
#endif
// #include "vector_vcl.h"

namespace imart
{

// Class image
template <typename type, typename container=vector_cpu<type>>
class image: public inherit_multiple<image<type,container>, space_object, data_object<type,container>>
{
public:
    // Type definitions
    using self    = image;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using space_object::dim;
    using space_object::size;
    using space_object::spacing;
    using space_object::direction;
    using data_object<type,container>::data;
    using data_object<type,container>::num_elements;
    using data_object<type,container>::get_data;

protected:
    // Type definitions
    // using iterator    = typename std::vector<type>::iterator;
    using container_pointer = typename data_object<type,container>::container_pointer;

    // ===========================================
    // Internal Variables
    // ===========================================
    int width;
    int height;
    int length;
    int channels;       // for now will only support one channel

    // Initialization methods as protected
    virtual void init(int w, int h, int channel);           // Initialize 2d
    virtual void init(int w, int h, int l, int channel);    // Initialize 3d
    virtual void copy_properties(const image & input); // copy only properties
    // virtual void allocate(int elements);

    // ===========================================
    // Interface Functions
    // ===========================================
    // Read methods protected
    void read_2d(std::string file_name);
    void read_3d(std::string file_name);
    // Write methods protected
    void write_2d(std::string file_name);
    void write_3d(std::string file_name);

public:

    // ===========================================
    // Constructor Functions
    // ===========================================
    // Constructors
    //! Constructor empty. Default a 2d image
    image();
    //! Constructor empty with dimension
    image(int d);
    //! Constructor using width and height.
    image(int w, int h);                           // Only 2d
    //! Constructor using width, height and length.
    image(int w, int h, int l);                    // Only 3d
    //! Constructor with vector size
    image(std::vector<int> shape);
    //! Constructor with vector size and channel
    image(std::vector<int> shape, int channel);
    //! Constructor with existing data.
    image(container_pointer buffer, int w, int h);        // Only 2d
    //! Constructor with existing data.
    image(container_pointer buffer, int w, int h, int l); // Only 3d
    //! Constructor with list
    image(std::initializer_list<type> list);
    //! Constructor with list
    image(std::vector<int> shape, std::initializer_list<type> list);

    //! Constructor copy
    image(const image & input);
    //! Constructor with file_name (call read()).
    // image(std::string file_name);
    //! Destructor
    // ~image();

    // ===========================================
    // Create Functions
    // ===========================================
    //! Full copy image
    virtual void clone_(const image & input);
    //! Duplicate data array for images
    virtual void copy_(const image & input);
    //! Imitate other image with a property copy
    virtual void mimic_(const image & input);
    // !Make the data equal. This function is to hold the memory pointer
    void equal(const image & input);
    //! Send data to cpu
    void to_cpu();

    // ===========================================
    // Get Functions
    // ===========================================
    //! Get the image width.
    int get_width() const;
    //! Get the image height.
    int get_height() const;
    //! Get the image length.
    int get_length() const;
    //! Get the image number of channels.
    int get_channels() const;
    //! Get the number of elements allocated.
    int get_total_elements() const;
    //! Get shared pointer count
    int get_ptr_count() const;
    //! Get shared pointer with image data
    container_pointer get_data() const;
    //! Get raw pointer to image data
    type * ptr() const;
    
    // //! Get the iterator begin
    // iterator begin() const;
    // //! Get the iterator end
    // iterator end() const;
    //! Assert function to process input image
    void assert_size(const image & input);

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);

    // ===========================================
    // Initialization Functions
    // ===========================================
    void zeros();
    void ones();
    void assign(type value);
    void random(float min=0.0, float max=1.0);
    // void identity(); //TODO

    // ===========================================
    // Overloading Operators
    // ===========================================
    // Access, does not support writing pixels
    type operator [] (int e);
    // type & operator () (int w, int h); //ONLY 2d***
    // type & operator () (int w, int h, int l); //ONLY 3d***
    
    // Equal
    image<type,container> & operator = (const image & input);
    
    // Image to Image
    image<type,container> operator + (const image & input);
    image<type,container> operator - (const image & input);
    image<type,container> operator * (const image & input);
    image<type,container> operator / (const image & input);
    image<type,container> operator ^ (const image & input);

    // Scalar
    image<type,container> operator + (type scalar);
    image<type,container> operator - (type scalar);
    image<type,container> operator * (type scalar);
    image<type,container> operator / (type scalar);
    image<type,container> operator ^ (type scalar);

    // Friend classes to support double side
    template<typename type_, typename container_>
    friend image<type_,container_> operator + (type_ scalar, image<type_,container_> & input);
    template<typename type_, typename container_>
    friend image<type_,container_> operator - (type_ scalar, const image<type_,container_> & input);
    template<typename type_, typename container_>
    friend image<type_,container_> operator * (type_ scalar, image<type_,container_> & input);
    template<typename type_, typename container_>
    friend image<type_,container_> operator / (type_ scalar, const image<type_,container_> & input);
    template<typename type_, typename container_>
    friend image<type_,container_> operator ^ (type_ scalar, const image<type_,container_> & input);
    
    // Conditions
    image<type,container> operator == (const image & input);
    image<type,container> operator > (const image & input);
    image<type,container> operator < (const image & input);
    image<type,container> operator >= (const image & input);
    image<type,container> operator <= (const image & input);

    image<type,container> operator == (type scalar);
    image<type,container> operator > (type scalar);
    image<type,container> operator < (type scalar);
    image<type,container> operator >= (type scalar);
    image<type,container> operator <= (type scalar);

    // replace indexes with new values
    void replace(const image<type,container> & idxs, const image & input);
    void replace(const image<type,container> & idxs, type value);

    // subregion
    pointer region(std::vector<int> offset, std::vector<int> out_size);

    // ===========================================
    // Reduction Functions
    // ===========================================
    int argmin();
    int argmax();
    type min();
    type max();
    type sum();
    // type prod();
    type dot(const image & input);

    type mean();
    type std_dev();

    // Matrix product //ONLY 2d***
    // image<type> _x_(const image<type> & input);

    // ===========================================
    // Interface Functions
    // ===========================================
    // interface with files using ITK
    void read(std::string file_name);
    void write(std::string file_name);

    //ITK interface, VTK?
    // template <type_itk>
    // void read_itk(itk::Image< type_itk, 2 > image_itk);
    // void write_itk();

    // TODO
    // create operator << to print info of image as image_info function [DONE] (inherited object)
    // class operations: +,-,*,/  [DONE]
    // initialize data with zeros, ones, random [DONE]
    // filters: normalize (0 to 1), padding, gaussian, convolution [not], gradient, fft [DONE]
    // scalar operations: scalar*Image [DONE]
    // functions in_place: transpose, add, substract, multiply, divide, pow
    // extra functions: copy, cast [DONE]
};

template<typename type>
using image_cpu = image<type,vector_cpu<type>>;

// template<typename type>
// using image_gpu = image<type,vector_opencl<type>>;

#ifdef IMART_WITH_OPENCL
template<typename type>
using image_opencl = image<type,vector_opencl<type>>;
#endif

#ifdef IMART_WITH_CUDA
template<typename type>
using image_cuda = image<type,vector_cuda<type>>;
#endif

// ===========================================
//      Functions of Class image
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename type, typename container>
image<type,container>::image()
{
    this->class_name = "image";
    init(0, 0, 1);
    // data.reset();
};

// Constructor
template <typename type, typename container>
image<type,container>::image(int d)
{
    // assert(d > 1 && d < 4);
    imart_assert(d > 1 && d < 4, "Dimension out of range");
    this->class_name = "image";
    if (d == 2) init(0, 0, 1);
    else if (d == 3) init(0, 0, 0, 1);
    // data.reset();
};

// Constructor
template <typename type, typename container>
image<type,container>::image(int w, int h)
{
    this->class_name = "image";
    init(w, h, 1);
};

// Constructor
template <typename type, typename container>
image<type,container>::image(int w, int h, int l)
{
    this->class_name = "image";
    init(w, h, l, 1);
};

// Constructor
template <typename type, typename container>
image<type,container>::image(std::vector<int> shape)
{
    this->class_name = "image";
    if (shape.size() == 2) init(shape[0], shape[1], 1);
    else if (shape.size() == 3) init(shape[0], shape[1], shape[2], 1);
    else ; // undefined for other dimensions
};

// Constructor
template <typename type, typename container>
image<type,container>::image(std::vector<int> shape, int channel)
{
    this->class_name = "image";
    if (shape.size() == 2) init(shape[0], shape[1], channel);
    else if (shape.size() == 3) init(shape[0], shape[1], shape[2], channel);
    else ; // undefined for other dimensions
};

// Constructor
template <typename type, typename container>
image<type,container>::image(container_pointer buffer, int w, int h)
{
    this->class_name = "image";
    // assert(buffer->size() == w*h);
    imart_assert(buffer->size() == w*h, "Buffer size different to width*height");
    init(w, h, 1);
    data.reset();
    data = buffer;
};

// Constructor
template <typename type, typename container>
image<type,container>::image(container_pointer buffer, int w, int h, int l)
{
    this->class_name = "image";
    // assert(buffer->size() == w*h*l);
    imart_assert(buffer->size() == w*h*l, "Buffer size different to width*height*length");
    init(w, h, l, 1);
    data.reset();
    data = buffer;
};

// // Constructor
template <typename type, typename container>
image<type,container>::image(std::initializer_list<type> list)
{
    this->class_name = "image";
    int s = list.size();
    init(s, 1, 1);
    data.reset();
    data = container::new_pointer(list);
};

template <typename type, typename container>
image<type,container>::image(std::vector<int> shape, std::initializer_list<type> list)
{
    this->class_name = "image";
    int d = shape.size();
    int s = list.size();

    int m = 1;
    for (int i=0; i<d; i++) m = m*shape[i];

    // printf("vector size %d\n",m);
    // printf("list size %d\n",s);
    imart_assert( m == s, "vector and initializer list have different size");
    
    if (d == 2)
        init(shape[0], shape[1], 1);
    else if (d == 3)
        init(shape[0], shape[1], shape[2], 1);
        
    data.reset();
    data = container::new_pointer(list);
};

// Constructor
template <typename type, typename container>
image<type,container>::image(const image & input)
{
    this->class_name = "image";
    clone_(input);
};

// Destructor
// template <typename type, typename container>
// image<type,container>::~image()
// {
//     ;//data.reset();
// };

// Empty init
template <typename type, typename container>
void image<type,container>::init(int w, int h, int channel)
{   
    channels = channel;
    width = w;
    height = h;
    length = 1;
    num_elements = width*height;

    space_object::init(2);
    data_object<type,container>::init(num_elements);   // allocation

    // this->size = std::vector<int>{width, height};
    space_object::set_size(std::vector<int>{width, height});
};

template <typename type, typename container>
void image<type,container>::init(int w, int h, int l, int channel)
{
    channels = channel;
    width = w;
    height = h;
    length = l;
    num_elements = width*height*length;

    space_object::init(3);                              // spatial properties
    data_object<type,container>::init(num_elements);    // allocate

    // this->size = std::vector<int>{width, height, length};
    space_object::set_size(std::vector<int>{width, height, length});
};

// Copy metadata
template <typename type, typename container>
void image<type,container>::copy_properties(const image & input)
{
    width = input.get_width();
    height = input.get_height();
    length = input.get_length();
    channels = input.get_channels();
    num_elements = input.get_total_elements();
};

// Full copy
template <typename type, typename container>
void image<type,container>::clone_(const image & input)
{
    copy_properties(input);
    space_object::clone_(input);
    data_object<type,container>::clone_(input);
};

// Point to the same data
template <typename type, typename container>
void image<type,container>::copy_(const image & input)
{
    copy_properties(input);
    space_object::copy_(input);
    data_object<type,container>::copy_(input);
};

template <typename type, typename container>
void image<type,container>::mimic_(const image & input)
{
    copy_properties(input);
    space_object::mimic_(input);
    data_object<type,container>::mimic_(input);
};

template <typename type, typename container>
void image<type,container>::equal(const image & input)
{
    assert_size(input);
    this->get_data()->equal(input.get_data());
};

template <typename type, typename container>
void image<type,container>::to_cpu()
{
    this->get_data()->to_cpu();
};

// template <typename type, typename container>
// typename image<type,container>::pointer image<type,container>::new_pointer()
// {
//     image<type,container>::pointer img = std::make_shared<image<type>>();
//     return img;
// };

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
int image<type,container>::get_width() const
{
    return width;
};

template <typename type, typename container>
int image<type,container>::get_height() const
{
    return height;
};

template <typename type, typename container>
int image<type,container>::get_length() const
{
    return length;
};

template <typename type, typename container>
int image<type,container>::get_channels() const
{
    return channels;
};

template <typename type, typename container>
int image<type,container>::get_total_elements() const
{
    return num_elements;
};

template <typename type, typename container>
int image<type,container>::get_ptr_count() const
{
    return data.use_count();
};

template <typename type, typename container>
typename image<type,container>::container_pointer image<type,container>::get_data() const
{
    return data;
};

template <typename type, typename container>
type * image<type,container>::ptr() const
{
    if (num_elements > 0)  { return data->data(); }
    else                   { return nullptr; };
    // return (*data).data();
};

// template <typename type, typename container>
// typename image<type,container>::iterator image<type,container>::begin() const
// {
//     return data->begin();
// };

// template <typename type, typename container>
// typename image<type,container>::iterator image<type,container>::end() const
// {
//     return data->end();
// };

// template <typename type, typename container>
// std::shared_ptr<std::vector<image<type>>> image<type,container>::get_grid()
// {
//     return grid;
// };

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string image<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Image Information";
    if (msg != "") { title = msg; };
    // Summary of the image information
    ss << object::info(title);
    ss << space_object::info(title);
    ss << "Pixel channels: \t" << channels << std::endl;
    ss << data_object<type,container>::info(title);
    return ss.str();
};

template <typename type, typename container>
std::string image<type,container>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; };
    // ss << data_object<type,container>::info_data("");
    std::vector<type> tmp = data->std_vector();
    type * p = tmp.data();

    if (this->dim == 2)
    {
        // std::cout << "Image data:" << std::endl;
        // std::cout << "["
        for(int i = 0; i < height; i++)
        {
            for(int j=0; j < width; j++)
            {
                ss << p[j+i*width] << " "; // valgrind error solved
            };
            ss << std::endl;
        };
        // std::cout << "]";
        ss << std::endl;
    };

    if (this->dim == 3)
    {
        // std::cout << "Image data:" << std::endl;
        // ss << "[";
        for(int k = 0; k < length; k++)
        {
            ss << "[ ";
            for(int i = 0; i < height; i++)
            {
                for(int j = 0; j < width; j++)
                {
                    ss << p[j + i*width + k*width*height] << " "; // valgrind error solved
                };
                if(i < height-1){ss << std::endl << "  ";};
            };
            ss << "]" << std::endl;
        };
        // ss << "]";
        ss << std::endl;
    };

    return ss.str();
};

// ===========================================
// Initialization Functions
// ===========================================
template <typename type, typename container>
void image<type,container>::zeros()
{
    data->zeros();
};

template <typename type, typename container>
void image<type,container>::ones()
{
    data->ones();
};

// Set all pixel to a fixed value
template <typename type, typename container>
void image<type,container>::assign(type value)
{
    data->assign(value);
};

template <typename type, typename container>
void image<type,container>::random(float min, float max)
{
    data->random();
};



// ===========================================
// Overloading Operators
// ===========================================
// Access
template <typename type, typename container>
type image<type,container>::operator [] (int e)
{
    imart_assert(e >= 0 && e < num_elements, "Image index out of bounds" );
    return data->operator[](e);
};

// template <typename type, typename container>
// type & image<type,container>::operator () (int w, int h)
// {
//     // assert(w < width and h < height);
//     return (*data)[w+h*width];
// };

// template <typename type, typename container>
// type & image<type,container>::operator () (int w, int h, int l)
// {
//     // assert(w < width and h < height and l < length);
//     return (*data)[w+h*width+l*height*width];
// };

template <typename type, typename container>
void image<type,container>::assert_size(const image & input)
{
    // assert(this->get_width()  == input.get_width());
    // assert(this->get_height() == input.get_height());
    // assert(this->get_length() == input.get_length());
    // assert(this->get_channels() == input.get_channels());
    imart_assert(this->get_width()  == input.get_width(), "Mismatch of image width");
    imart_assert(this->get_height() == input.get_height(), "Mismatch of image height");
    imart_assert(this->get_length() == input.get_length(), "Mismatch of image length");
    imart_assert(this->get_channels() == input.get_channels(), "Mismatch of image channels");
    return;
};

// Equal
template <typename type, typename container>
image<type,container> & image<type,container>::operator = (const image & input)
{
    // delete &data;
    copy_(input);
    return *this;
};

// Image to Image
template <typename type, typename container>
image<type,container> image<type,container>::operator + (const image & input)
{
    assert_size(input);
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) + *(input.get_data()));
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator - (const image & input)
{
    assert_size(input);
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) - *(input.get_data()));
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator * (const image & input)
{
    assert_size(input);
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) * *(input.get_data()));
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator / (const image & input)
{
    assert_size(input);
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) / *(input.get_data()));
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator ^ (const image & input)
{
    assert_size(input);
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) ^ *(input.get_data()));
    return *output;
};

// Scalar right hand side
template <typename type, typename container>
image<type,container> image<type,container>::operator + (type scalar)
{
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) + scalar);
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator - (type scalar)
{
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) - scalar);
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator * (type scalar)
{
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) * scalar);
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator / (type scalar)
{
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) / scalar);
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator ^ (type scalar)
{
    image<type,container>::pointer output = this->mimic();
    output->set_data(*(this->get_data()) ^ scalar);
    return *output;
};

// Scalar left hand side defined in header due to use of friend functions
template <typename type, typename container>
image<type,container> operator + (type scalar, image<type,container> & input)
{
    return input + scalar;
};

template <typename type, typename container>
image<type,container> operator - (type scalar, const image<type,container> & input)
{
    typename image<type,container>::pointer output = input.mimic();
    output->set_data(scalar - *(input.get_data()));
    return *output;
};

template <typename type, typename container>
image<type,container> operator * (type scalar, image<type,container> & input)
{
    return input * scalar;
};

template <typename type, typename container>
image<type,container> operator / (type scalar, const image<type,container> & input)
{
    typename image<type,container>::pointer output = input.mimic();
    output->set_data(scalar / *(input.get_data()));
    return *output;
};

template <typename type, typename container>
image<type,container> operator ^ (type scalar, const image<type,container> & input)
{
    typename image<type,container>::pointer output = input.mimic();
    output->set_data(scalar ^ *(input.get_data()));
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator == (const image<type,container> & input)
{
    auto output = image<type,container>::new_pointer(input.get_size());
    output->set_data( *(this->get_data()) == *(input.get_data()) );
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator > (const image<type,container> & input)
{
    auto output = image<type,container>::new_pointer(input.get_size());
    output->set_data( *(this->get_data()) > *(input.get_data()) );
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator < (const image<type,container> & input)
{
    auto output = image<type,container>::new_pointer(input.get_size());
    output->set_data( *(this->get_data()) < *(input.get_data()) );
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator >= (const image<type,container> & input)
{
    auto output = image<type,container>::new_pointer(input.get_size());
    output->set_data( *(this->get_data()) >= *(input.get_data()) );
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator <= (const image<type,container> & input)
{
    auto output = image<type,container>::new_pointer(input.get_size());
    output->set_data( *(this->get_data()) <= *(input.get_data()) );
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator == (type scalar)
{
    auto output = image<type,container>::new_pointer(this->get_size());
    output->set_data(*(this->get_data()) == scalar);
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator > (type scalar)
{
    auto output = image<type,container>::new_pointer(this->get_size());
    output->set_data(*(this->get_data()) > scalar);
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator < (type scalar)
{
    auto output = image<type,container>::new_pointer(this->get_size());
    output->set_data(*(this->get_data()) < scalar);
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator >= (type scalar)
{
    auto output = image<type,container>::new_pointer(this->get_size());
    output->set_data(*(this->get_data()) >= scalar);
    return *output;
};

template <typename type, typename container>
image<type,container> image<type,container>::operator <= (type scalar)
{
    auto output = image<type,container>::new_pointer(this->get_size());
    output->set_data(*(this->get_data()) <= scalar);
    return *output;
};

template <typename type, typename container>
void image<type,container>::replace(const image<type,container> & idxs, const image & input)
{
    this->get_data()->replace( *(idxs.get_data()), *(input.get_data()) );
    return;
};

template <typename type, typename container>
void image<type,container>::replace(const image<type,container> & idxs, type value)
{
    this->get_data()->replace( *(idxs.get_data()), value );
    return;
};

template <typename type, typename container>
typename image<type,container>::pointer image<type,container>::region(std::vector<int> offset, std::vector<int> out_size)
{
    auto output = image<type,container>::new_pointer(out_size);
    std::vector<int> ref_size = this->get_size();
    output->set_data( this->get_data()->region(ref_size, offset, out_size) );
    return output;
};


// ===========================================
// Functions
// ===========================================
// Reductions operations
// Reductions not using OpenMP. They are slower.
template <typename type, typename container>
type image<type,container>::min()
{
    return data->min();
};

template <typename type, typename container>
type image<type,container>::max()
{
    return data->max();
};

template <typename type, typename container>
int image<type,container>::argmin()
{
    return data->argmin();
};

template <typename type, typename container>
int image<type,container>::argmax()
{
    return data->argmax();
};

template <typename type, typename container>
type image<type,container>::sum()
{
    return data->sum();
};

// template <typename type, typename container>
// type image<type,container>::prod()
// {
    
// };

// Vectorial dot product. Verify the same number of elements, then product and reduce
template <typename type, typename container>
type image<type,container>::dot(const image & input)
{
    assert_size(input);
    return data->dot(*(input.get_data()));
};

template <typename type, typename container>
type image<type,container>::mean()
{
    return type(1.0/get_total_elements())*data->sum();
};

template <typename type, typename container>
type image<type,container>::std_dev()
{   
    type mu = mean();
    auto sub = *(data) - mu;
    sub = (*sub)*(*sub);
    return std::sqrt( sub->sum()/ (type(get_total_elements()-1)) );
};

// Matrix product
// template <typename type, typename container>
// image<type> image<type,container>::_x_(const image<type> & input)
// {
//     using Map = Eigen::Map<Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
//     image<type> result(input.width, height);
//     // std::shared_ptr<image<type>> result = 
//     // std::shared_ptr<image<type>>(new image<type>(height, input.width));
//     // result.print();
//     // std::cout << "ptr: " << result.get_data() << std::endl;

//     Map A(ptr(), height, width);
//     Map B(input.ptr(), input.get_height(), input.get_width());
//     Map C(result.ptr(), result.get_height(), result.get_width());

//     // std::cout << A << std::endl;
//     // std::cout << std::endl;
//     // std::cout << B << std::endl;
//     // std::cout << std::endl;
//     // std::cout << A*B << std::endl;
//     // std::cout << std::endl;

//     C = (A*B);
//     // C = (A.transpose()*B.transpose()).transpose(); // transpose function doesn't cost anything (average=12ns)
//     return result;
// };

// ===========================================
// Interface Functions
// ===========================================
template <typename type, typename container>
void image<type,container>::read(std::string file_name)
{
    if (this->get_dimension() == 2) { read_2d(file_name); }
    else if (this->get_dimension() == 3) { read_3d(file_name); };
};

template <typename type, typename container>
void image<type,container>::read_2d(std::string file_name)
{
    const int d = 2;
    using itkImageType = itk::Image<type, d>;
    using itkReaderType = itk::ImageFileReader<itkImageType>;
    
    // Initialize itk objects
    typename itkImageType::Pointer image_itk = itkImageType::New();
    typename itkReaderType::Pointer reader = itkReaderType::New();

    // Set the image filename itk
    reader->SetFileName(file_name);
    reader->Update();

    // Read the image from reader
    image_itk = reader->GetOutput();
    // std::cout << image_itk;

    typename itkImageType::RegionType region = image_itk->GetLargestPossibleRegion();
    typename itkImageType::SizeType itksize = region.GetSize();
    type * p = image_itk->GetBufferPointer();
    int w = itksize[0];
    int h = itksize[1];

    init(w, h, 1);
    // data.reset();
    // allocate(num_elements);

    // Copy of the image
    data->read_ram(p, num_elements);
    // #pragma omp parallel for // ** Better times without omp
    // for(int k=0; k<num_elements; k++)
    // {
    //     (*data)[k] = *(p+k);
    // };

    typename itkImageType::SpacingType spacing_itk = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk = image_itk->GetDirection();

    std::vector<double> spacing_im(d,0), origin_im(d,0), direction_im(d*d,0);
    for (int i = 0; i < d; i++)
    {
        spacing_im[i] = spacing_itk[i];
        origin_im[i] = origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction_im[c] = direction_itk[i][j];
            c++;
        };
    };

    this->set_spacing(spacing_im);
    this->set_origin(origin_im);
    this->set_direction(direction_im);
};

template <typename type, typename container>
void image<type,container>::read_3d(std::string file_name)
{
    const int d = 3;
    using itkImageType = itk::Image<type, d>;
    using itkReaderType = itk::ImageFileReader<itkImageType>;
    
    // Initialize itk objects
    typename itkImageType::Pointer image_itk = itkImageType::New();
    typename itkReaderType::Pointer reader = itkReaderType::New();

    // Set the image filename itk
    reader->SetFileName(file_name);
    reader->Update();

    // Read the image from reader
    image_itk = reader->GetOutput();
    // std::cout << image_itk;

    typename itkImageType::RegionType region = image_itk->GetLargestPossibleRegion();
    typename itkImageType::SizeType itksize = region.GetSize();
    type * p = image_itk->GetBufferPointer();
    int w = itksize[0];
    int h = itksize[1];
    int l = itksize[2];

    init(w, h, l, 1);
    // data.reset();
    // allocate(num_elements);

    // Copy of the image
    data->read_ram(p, num_elements);
    // #pragma omp parallel for // ** Better times without omp
    // for(int k=0; k<num_elements; k++)
    // {
    //     (*data)[k] = *(p+k);
    // };

    typename itkImageType::SpacingType spacing_itk = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk = image_itk->GetDirection();

    std::vector<double> spacing_im(d,0), origin_im(d,0), direction_im(d*d,0);
    for (int i = 0; i < d; i++)
    {
        spacing_im[i] = spacing_itk[i];
        origin_im[i] = origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction_im[c] = direction_itk[i][j];
            c++;
        };
    };

    this->set_spacing(spacing_im);
    this->set_origin(origin_im);
    this->set_direction(direction_im);
};

// template <size_t type_itk>
// void image::read_itk(itk::Image< type_itk, 2 > image_itk)
// {
//     template <size_t type_itk> 
//     using ImageType = itk::Image< type_itk, 2 >;
//     // typedef itk::Image< typename type_itk, 2 > ImageType;

//     ImageType::RegionType region = image_itk->GetLargestPossibleRegion();
//     ImageType::SizeType size = region.GetSize();
//     type * p = image_itk->GetBufferPointer();
//     int w = size[0];
//     int h = size[1];

//     image(w,h); // empty current image and create new

//     // Copy of the image
//     for(int k=0; k<w*h; k++)
//     {
//         *(data+k) = *(p+k);
//     };
// };

template <typename type, typename container>
void image<type,container>::write(std::string file_name)
{
    if (this->get_dimension() == 2) { write_2d(file_name); }
    else if (this->get_dimension() == 3) { write_3d(file_name); };
};

template <typename type, typename container>
void image<type,container>::write_2d(std::string file_name)
{
    const int d = 2;
    using itkImageType = itk::Image<type, d>;
    using itkWriterType = itk::ImageFileWriter<itkImageType>;

    typename itkImageType::Pointer image_itk = itkImageType::New();
    typename itkWriterType::Pointer writer = itkWriterType::New();

    typename itkImageType::RegionType region;
    typename itkImageType::IndexType  start;
    typename itkImageType::SizeType itksize;
    
    start[0] = 0;
    start[1] = 0;
    itksize[0] = get_width();
    itksize[1] = get_height();

    region.SetSize(itksize);
    region.SetIndex(start);

    image_itk->SetRegions(region);
    image_itk->Allocate();

    type * p = image_itk->GetBufferPointer();

    data->write_ram(p,num_elements);
    // for(int k=0; k<num_elements; k++)
    // {
    //     *(p+k) = (*data)[k];
    // };
    
    // Writing metadata
    std::vector<double> spacing_im = this->get_spacing();
    std::vector<double> origin_im = this->get_origin();
    std::vector<double> direction_im = this->get_direction();
    
    typename itkImageType::SpacingType spacing_itk;// = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk;// = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk;// = image_itk->GetDirection();
    
    for (int i = 0; i < d; i++)
    {
        spacing_itk[i] = spacing_im[i];
        origin_itk[i] = origin_im[i];
        // std::cout << origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction_itk[i][j] = direction_im[c];
            c++;
        };
    };

    image_itk->SetSpacing(spacing_itk);
    image_itk->SetOrigin(origin_itk);
    image_itk->SetDirection(direction_itk);

    writer->SetFileName(file_name);
    writer->SetInput(image_itk);

    try
    {
        writer->Update();
    }
    catch (itk::ExceptionObject & error)
    {
        std::cerr << "Error: " << error << std::endl;
    }
};

template <typename type, typename container>
void image<type,container>::write_3d(std::string file_name)
{
    const int d = 3;
    using itkImageType = itk::Image<type, d>;
    using itkWriterType = itk::ImageFileWriter<itkImageType>;

    typename itkImageType::Pointer image_itk = itkImageType::New();
    typename itkWriterType::Pointer writer = itkWriterType::New();

    typename itkImageType::RegionType region;
    typename itkImageType::IndexType  start;
    typename itkImageType::SizeType itksize;
    
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    itksize[0] = get_width();
    itksize[1] = get_height();
    itksize[2] = get_length();

    region.SetSize(itksize);
    region.SetIndex(start);

    image_itk->SetRegions(region);
    image_itk->Allocate();

    type * p = image_itk->GetBufferPointer();

    data->write_ram(p,num_elements);
    // for(int k=0; k<num_elements; k++)
    // {
    //     *(p+k) = (*data)[k];
    // };
    
    // Writing metadata
    std::vector<double> spacing_im = this->get_spacing();
    std::vector<double> origin_im = this->get_origin();
    std::vector<double> direction_im = this->get_direction();
    
    // When I uncomment this lines spacing in z goes to 10^14???

    typename itkImageType::SpacingType spacing_itk;// = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk;// = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk;// = image_itk->GetDirection();

    for (int i = 0; i < d; i++)
    {
        spacing_itk[i] = spacing_im[i];
        origin_itk[i] = origin_im[i];
        // std::cout << origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction_itk[i][j] = direction_im[c];
            c++;
        };
    };

    image_itk->SetSpacing(spacing_itk);
    image_itk->SetOrigin(origin_itk);
    image_itk->SetDirection(direction_itk);
    
    writer->SetFileName(file_name);
    writer->SetInput(image_itk);

    try
    {
        writer->Update();
    }
    catch (itk::ExceptionObject & error)
    {
        std::cerr << "Error: " << error << std::endl;
    }
};

}; //end namespace

// Template constructions
// template class image<unsigned char>;  // 1 byte
// template class image<unsigned short>; // 2 byte
// template class image<short>;          // 2 byte
// template class image<unsigned int>;   // 4 byte
// template class image<int>;            // 4 byte
// template class image<float>;          // 4 byte
// template class image<double>;         // 8 byte

#endif