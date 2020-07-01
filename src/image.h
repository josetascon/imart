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
#include <random>       // random
#include <cassert>      // assert
#include <cmath>        // math (pow)
#include <algorithm>    // std::max std::min
#include <complex>      // std::complex

#include <fftw3.h>      // fft library

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
#include "vector_ocl.h"
#include "vector_vcl.h"

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

protected:
    // Type definitions
    // using iterator    = typename std::vector<type>::iterator;
    using pixels4 = std::unique_ptr<std::array<type,4>>;
    using pixels8 = std::unique_ptr<std::array<type,8>>;
    using container_pointer = typename data_object<type,container>::container_pointer;

    // ===========================================
    // Internal Variables
    // ===========================================
    int width;
    int height;
    int length;
    int channels;       // for now will only support one channel

    // std::vector<int> size;
    // std::vector<type> spacing;
    // std::vector<type> origin;
    // std::vector<type> direction; //not supported with interpolation

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
    image(std::vector<int> size, int channel);
    // //! Constructor with existing data.
    image(container_pointer buffer, int w, int h);        // Only 2d
    // //! Constructor with existing data.
    image(container_pointer buffer, int w, int h, int l); // Only 3d
    // //! Constructor with list
    image(std::initializer_list<type> list);

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

    // typename image<type,container>::pointer new_pointer();

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
    void identity(); //TODO

    // ===========================================
    // Overloading Operators
    // ===========================================
    // Access, does not support writing pixels
    type & operator [] (int e);
    // type & operator () (int w, int h); //ONLY 2d***
    // type & operator () (int w, int h, int l); //ONLY 3d***

    // Access
    // pixels2 neighbors(int e);
    pixels4 neighbors4(int e);
    pixels8 neighbors8(int e);
    
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
    
    // ===========================================
    // Functions
    // ===========================================
    // Matrix product //ONLY 2d***
    // image<type> _x_(const image<type> & input);
    
    // reduction function
    type min();
    type max();
    type sum();
    // type prod();
    type dot(const image & input);

    template<typename type_, typename container_>
    friend image<type_,container_> normalize(const image<type_,container_> & input, type_ min, type_ max);

    template<typename type_out, typename container_, typename type_in>
    friend image<type_out,container_> cast(const image<type_in,container_> & input);

    // image utils functions
    template<typename type_, typename container_>
    friend image<type_,container_> pad(const image<type_,container_> & input, std::vector<int> pre, std::vector<int> post);

    /*
    template<typename pixel_t>
    friend image<pixel_t> unpad(const image<pixel_t> & input, std::vector<int> pre, std::vector<int> post);
    
    // friend functions
    // template<typename pixel_t>
    // friend image<std::complex<pixel_t>> real2complex(const image<pixel_t> & input);
    template<typename pixel_t>
    friend image<pixel_t> complex2real(const image<std::complex<pixel_t>> & input);
    
    // fft in 2d only now
    template<typename pixel_t>
    friend image<std::complex<pixel_t>> fft(const image<pixel_t> & input);
    template<typename pixel_t>
    friend image<std::complex<pixel_t>> ifft(const image<std::complex<pixel_t>> & input);
    
    template<typename pixel_t>
    friend typename image<pixel_t>::vector gradient(const image<pixel_t> & input);
    // std::vector<image<type>> gradient();
    */
    // ===========================================
    // Interface Functions
    // ===========================================
    // interface with files using ITK
    void read(std::string file_name);
    void write(std::string file_name);

    //ITK interface, VTK? eigen?
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

template<typename type>
using image_gpu = image<type,vector_ocl<type>>;


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
    data.reset();
};

// Constructor
template <typename type, typename container>
image<type,container>::image(int d)
{
    assert(d > 1 && d < 4);
    this->class_name = "image";
    if (d == 2) init(0, 0, 1);
    if (d == 3) init(0, 0, 0, 1);
    data.reset();
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
image<type,container>::image(std::vector<int> size, int channel)
{
    this->class_name = "image";
    if (size.size() == 2) init(size[0], size[1], channel);
    else if (size.size() == 3) init(size[0], size[1], size[2], channel);
    else ; // undefined for other dimensions
};

// Constructor
template <typename type, typename container>
image<type,container>::image(container_pointer buffer, int w, int h)
{
    this->class_name = "image";
    assert(buffer->size() == w*h);
    init(w, h, 1);
    data.reset();
    data = buffer;
};

// Constructor
template <typename type, typename container>
image<type,container>::image(container_pointer buffer, int w, int h, int l)
{
    this->class_name = "image";
    assert(buffer->size() == w*h*l);
    init(w, h, l, 1);
    data.reset();
    data = buffer;
};

// // Constructor
template <typename type, typename container>
image<type,container>::image(std::initializer_list<type> list)
{
    int s = list.size();
    init(s, 1, 1);
    data.reset();
    data = container::new_pointer(list);
}

// Constructor
template <typename type, typename container>
image<type,container>::image(const image & input)
{
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

    this->size = std::vector<int>{width, height};
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

    this->size = std::vector<int>{width, height, length};
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
type & image<type,container>::operator [] (int e)
{
    // assert(e < num_elements);
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
    assert(this->get_width()  == input.get_width());
    assert(this->get_height() == input.get_height());
    assert(this->get_length() == input.get_length());
    assert(this->get_channels() == input.get_channels());
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
    output->set_data(scalar - *(input.get_data()));
    return *output;
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

template<typename type, typename container>
image<type,container> normalize(const image<type,container> & input, type min = 0.0, type max = 1.0)
{
    typename image<type,container>::pointer output = input.mimic();
    output->set_data((input.get_data()->normalize(min,max)));
    return *output;
};

template<typename type_out, typename container_out, typename type_in, typename container_in>
image<type_out,container_out> cast(const image<type_in,container_in> & input)
{
    typename image<type_out,container_out>::pointer output = image<type_out,container_out>::new_pointer(input.get_size(),1);
    output->set_data(input.get_data()->template cast<type_out>());
    output->set_sod_parameters(input.get_spacing(), input.get_origin(), input.get_direction());
    return *output;
};


template<typename type, typename container>
image<type,container> pad(const image<type,container> & input, std::vector<int> pre, std::vector<int> post)
{
    int w = input.get_width();
    int h = input.get_height();
    int l = input.get_length();
    std::vector<int> extra(pre.size(),0);
    for (int i = 0; i < extra.size(); i++){ extra[i] = pre[i]+post[i]; };
    
    typename image<type,container>::pointer output;
    if (input.get_dimension() == 2){ output = image<type,container>::new_pointer(w+extra[0], h+extra[1]); };
    if (input.get_dimension() == 3){ output = image<type,container>::new_pointer(w+extra[0], h+extra[1], l+extra[2]); };

    container::pad(input.get_data(), output->get_data(), input.get_size(), pre, post);
    output->set_sod_parameters(input.get_spacing(), input.get_origin(), input.get_direction());
    return *output;
};

template<typename type, typename container>
image<type,container> unpad(const image<type,container> & input, std::vector<int> pre, std::vector<int> post)
{
    int w = input.get_width();
    int h = input.get_height();
    int l = input.get_length();
    std::vector<int> extra(pre.size(),0);
    for (int i = 0; i < extra.size(); i++){ extra[i] = pre[i]+post[i]; };
    
    typename image<type,container>::pointer output;
    if (input.get_dimension() == 2){ output = image<type,container>::new_pointer(w-extra[0], h-extra[1]); };
    if (input.get_dimension() == 3){ output = image<type,container>::new_pointer(w-extra[0], h-extra[1], l-extra[2]); };
    
    container::unpad(input.get_data(), output->get_data(), output->get_size(), pre, post);
    output->set_sod_parameters(input.get_spacing(), input.get_origin(), input.get_direction());
    return *output;
};

/*
template<typename type>
image<type> complex2real(const image<std::complex<type>> & input)
{
    int d = input.get_dimension();
    int N = input.get_total_elements();
    std::complex<type> * p1 = input.ptr();

    image<type> result(input.get_size());
    
    type * p2 = result.ptr();

    for (int i = 0; i < N; i++)
    {
        p2[i] = std::real(p1[i]);
    };
    return result;
};

template <typename type, typename container>
image<std::complex<type>> fft(const image<type> & input)
{
    int w = input.get_width();
    int h = input.get_height();
    int N = input.get_total_elements();
    fftw_complex in[N], out[N];
    type * p1 = input.ptr();

    for (int i = 0; i < N; i++)
    {
        in[i][0] = p1[i];
        in[i][1] = 0;
    };

    fftw_plan p_fft;
    p_fft = fftw_plan_dft_2d(h, w, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p_fft);
    fftw_destroy_plan(p_fft);

    image<std::complex<type>> result(w,h);
    std::complex<type> * p2 = result.ptr();    
    for (int i = 0; i < N; i++)
    {
        p2[i] = std::complex(out[i][0],out[i][1]);
    };

    return result;
};

template <typename type, typename container>
image<std::complex<type>> ifft(const image<std::complex<type>> & input)
{
    int w = input.get_width();
    int h = input.get_height();
    int N = input.get_total_elements();
    fftw_complex in[N], out[N];
    std::complex<type> * p1 = input.ptr();

    std::complex<type> v;
    for (int i = 0; i < N; i++)
    {
        v = p1[i];
        in[i][0] = v.real();
        in[i][1] = v.imag();
    };

    fftw_plan p_fft;
    p_fft = fftw_plan_dft_2d(h, w, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p_fft);
    fftw_destroy_plan(p_fft);

    image<std::complex<type>> result(w,h);
    std::complex<type> * p2 = result.ptr();    
    for (int i = 0; i < N; i++)
    {
        p2[i] = std::complex(out[i][0]/N,out[i][1]/N);
    };

    return result;
};

template<typename type>
typename image<type,container>::vector gradient(const image<type> & input)
{
    std::vector<int> none{0,0};    
    std::vector<int> extra{1,1};
    image<type> input_pad = pad(input, none, extra);
    // input_pad.print_data();
    int w = input_pad.get_width();
    int h = input_pad.get_height();

    // derivative x
    image<type> dx(w,h);
    type * pdx = dx.ptr();
    pdx[0] = 0.5;
    pdx[1] = 0;
    pdx[2] = -0.5;
    // dx.print_data();

    // derivative y
    image<type> dy(w,h);
    type * pdy = dy.ptr();
    pdy[0] = 0.5;
    pdy[w] = 0.0;
    pdy[2*w] = -0.5;
    // dy.print_data();

    typename image<type,container>::vector grad(2);
    typename image<type,container>::pointer didx(new image<type>(w,h));
    typename image<type,container>::pointer didy(new image<type>(w,h));
    grad[0] = didx;
    grad[1] = didy;

    // image<type> didx(w,h);
    image<type> didx_ = complex2real(ifft(fft(input_pad)*fft(dx)));
    *grad[0] = unpad(didx_, std::vector<int>{1,0}, std::vector<int>{0,1});
    // image<type> didx = unpad(didx_, std::vector<int>{1,0}, std::vector<int>{0,1});

    // image<type> didy(w,h);
    image<type> didy_ = complex2real(ifft(fft(input_pad)*fft(dy)));
    *grad[1] = unpad(didy_, std::vector<int>{0,1}, std::vector<int>{1,0});
    // image<type> didy = unpad(didy_, std::vector<int>{0,1}, std::vector<int>{1,0});

    // return didx;
    // return didy;
    return grad;
};
*/

// ===========================================
// Interface Functions
// ===========================================
template <typename type, typename container>
void image<type,container>::read(std::string file_name)
{
    if (this->get_dimension() == 2) { read_2d(file_name); };
    if (this->get_dimension() == 3) { read_3d(file_name); };
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

    std::vector<double> spacing(d,0), origin(d,0), direction(d*d,0);
    for (int i = 0; i < d; i++)
    {
        spacing[i] = spacing_itk[i];
        origin[i] = origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction[c] = direction_itk[i][j];
            c++;
        };
    };

    this->set_spacing(spacing);
    this->set_origin(origin);
    this->set_direction(direction);
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

    std::vector<double> spacing(d,0), origin(d,0), direction(d*d,0);
    for (int i = 0; i < d; i++)
    {
        spacing[i] = spacing_itk[i];
        origin[i] = origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction[c] = direction_itk[i][j];
            c++;
        };
    };

    this->set_spacing(spacing);
    this->set_origin(origin);
    this->set_direction(direction);
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
//     int h = size[0];

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
    if (this->get_dimension() == 2) { write_2d(file_name); };
    if (this->get_dimension() == 3) { write_3d(file_name); };
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
    typename itkImageType::SizeType size;
    
    start[0] = 0;
    start[1] = 0;
    size[0] = get_width();
    size[1] = get_height();

    region.SetSize(size);
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
    std::vector<double> spacing = this->get_spacing();
    std::vector<double> origin = this->get_origin();
    std::vector<double> direction = this->get_direction();
    
    typename itkImageType::SpacingType spacing_itk;// = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk;// = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk;// = image_itk->GetDirection();
    
    for (int i = 0; i < d; i++)
    {
        spacing_itk[i] = spacing[i];
        origin_itk[i] = origin[i];
        // std::cout << origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction_itk[i][j] = direction[c];
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
    typename itkImageType::SizeType size;
    
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    size[0] = get_width();
    size[1] = get_height();
    size[2] = get_length();

    region.SetSize(size);
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
    std::vector<double> spacing = this->get_spacing();
    std::vector<double> origin = this->get_origin();
    std::vector<double> direction = this->get_direction();
    
    // When I uncomment this lines spacing in z goes to 10^14???

    typename itkImageType::SpacingType spacing_itk;// = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk;// = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk;// = image_itk->GetDirection();

    for (int i = 0; i < d; i++)
    {
        spacing_itk[i] = spacing[i];
        origin_itk[i] = origin[i];
        // std::cout << origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction_itk[i][j] = direction[c];
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
typename image<type,container>::pixels4 image<type,container>::neighbors4(int e)
{
    auto vv = std::make_unique<std::array<type,4>>();
    // type * p = this->ptr();
    int w = this->width;
    (*vv)[0] = data->operator[](e);
    (*vv)[1] = data->operator[](e+1);
    (*vv)[2] = data->operator[](e+w);
    (*vv)[3] = data->operator[](e+w+1);
    // std::cout << "px:" << p[e] << " ";
    return vv;
};

template <typename type, typename container>
typename image<type,container>::pixels8 image<type,container>::neighbors8(int e)
{
    auto vv = std::make_unique<std::array<type,8>>();
    // type * p = this->ptr();
    int w = this->width;
    int h = this->height;

    (*vv)[0] = data->operator[](e);
    (*vv)[1] = data->operator[](e+1);
    (*vv)[2] = data->operator[](e+w);
    (*vv)[3] = data->operator[](e+w+1);
    
    (*vv)[4] = data->operator[](e+w*h);
    (*vv)[5] = data->operator[](e+1+w*h);
    (*vv)[6] = data->operator[](e+w+w*h);
    (*vv)[7] = data->operator[](e+w+1+w*h);
    // std::cout << "px:" << p[e] << " ";
    return vv;
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