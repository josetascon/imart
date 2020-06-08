/*
* @Author: jose
* @Date:   2019-11-05 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-14 17:39:46
*/

#ifndef __IMAGE_BASE_H__
#define __IMAGE_BASE_H__

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
#include <eigen3/Eigen/Core>

// local libs
#include "space_object.h"

namespace imart
{

// Class image_base
template <typename pixel_type>
class image_base: public space_object<pixel_type>
{
public:
    // Type definitions
    using self    = image_base;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using space_object<pixel_type>::dim;
    using space_object<pixel_type>::size;
    using space_object<pixel_type>::spacing;
    using space_object<pixel_type>::direction;

protected:
    // Type definitions
    using iterator    = typename std::vector<pixel_type>::iterator;
    using ptr_vector  = std::shared_ptr<std::vector<pixel_type>>;
    using ptr_pixels4 = std::unique_ptr<std::array<pixel_type,4>>;
    using ptr_pixels8 = std::unique_ptr<std::array<pixel_type,8>>;

    // ===========================================
    // Internal Variables
    // ===========================================
    int width;
    int height;
    int length;
    int channels;       // for now will only support one channel
    int num_elements;

    // std::vector<int> size;
    // std::vector<pixel_type> spacing;
    // std::vector<pixel_type> origin;
    // std::vector<pixel_type> direction; //not supported with interpolation

    // Image data
    ptr_vector data;

    // Initialization methods as protected
    virtual void init(int w, int h);            // Only 2d
    virtual void init(int w, int h, int l);     // Only 3d
    virtual void copy_properties(const image_base<pixel_type> & input); // copy only properties
    virtual void allocate(int elements);

    // Read methods protected
    void read_2d(std::string file_name);
    void read_3d(std::string file_name);
    // Write methods protected
    void write_2d(std::string file_name);
    void write_3d(std::string file_name);

public:

    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    //! Constructor empty. Default a 2d image
    image_base();
    //! Constructor empty with dimension
    image_base(int d);
    //! Constructor using width and height.
    image_base(int w, int h);                           // Only 2d
    //! Constructor using width, height and length.
    image_base(int w, int h, int l);                    // Only 3d
    //! Constructor with vector size
    image_base(std::vector<int> size);
    //! Constructor with existing data.
    image_base(ptr_vector buffer, int w, int h);        // Only 2d
    //! Constructor with existing data.
    image_base(ptr_vector buffer, int w, int h, int l); // Only 3d
    //! Constructor with list
    image_base(std::initializer_list<pixel_type> list);

    //! Constructor copy
    image_base(const image_base<pixel_type> & input);
    //! Constructor with file_name (call read()).
    image_base(std::string file_name);
    //! Destructor
    // ~image_base();

    //! New pointer
    // static image_base<pixel_type>::pointer new_pointer();
    template<typename... ARGS>
    static pointer new_pointer(const ARGS&... args);
    // { return std::make_shared<image_base>(args...);};

    //! Full copy image
    virtual void clone(const image_base & input);
    //! Duplicate data array for images
    virtual void copy(const image_base & input);
    //! Imitate other image with a property copy
    virtual void mimic(const image_base & input);

    // typename image_base<pixel_type>::pointer new_pointer();

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
    ptr_vector get_data() const;
    //! Get raw pointer to image data
    pixel_type * ptr() const;

    iterator begin() const;

    iterator end() const;

    void assert_size(const image_base<pixel_type> & input);

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg); //ONLY 2d***
    void print_ptr_count();

    // ===========================================
    // Initialization Functions
    // ===========================================
    void zeros();
    void ones();
    void identity(); //TODO
    void random(float min=0.0, float max=1.0);
    void fill(pixel_type value);

    // ===========================================
    // Overloading Operators
    // ===========================================
    // Access, does not support writing pixels
    pixel_type & operator () (int e);
    pixel_type & operator () (int w, int h); //ONLY 2d***
    pixel_type & operator () (int w, int h, int l); //ONLY 3d***

    // Equal
    image_base<pixel_type> & operator = (const image_base<pixel_type> & input);
    
    // Image to Image
    image_base<pixel_type> operator + (const image_base<pixel_type> & input);
    image_base<pixel_type> operator - (const image_base<pixel_type> & input);
    image_base<pixel_type> operator * (const image_base<pixel_type> & input);
    image_base<pixel_type> operator / (const image_base<pixel_type> & input);
    image_base<pixel_type> operator ^ (const image_base<pixel_type> & input);

    // Scalar
    image_base<pixel_type> operator + (pixel_type scalar);
    image_base<pixel_type> operator - (pixel_type scalar);
    image_base<pixel_type> operator * (pixel_type scalar);
    image_base<pixel_type> operator / (pixel_type scalar);
    image_base<pixel_type> operator ^ (pixel_type scalar);

    // Friend classes to support double side
    template<typename pixel_t>
    friend image_base<pixel_t> operator + (pixel_t scalar, image_base<pixel_t> & input);
    template<typename pixel_t>
    friend image_base<pixel_t> operator - (pixel_t scalar, const image_base<pixel_t> & input);
    template<typename pixel_t>
    friend image_base<pixel_t> operator * (pixel_t scalar, image_base<pixel_t> & input);
    template<typename pixel_t>
    friend image_base<pixel_t> operator / (pixel_t scalar, const image_base<pixel_t> & input);

    // ===========================================
    // Functions
    // ===========================================
    // Matrix product //ONLY 2d***
    image_base<pixel_type> _x_(const image_base<pixel_type> & input);
    
    // reduction function
    pixel_type min();
    pixel_type max();
    pixel_type sum();
    pixel_type prod();
    pixel_type dot(const image_base<pixel_type> & input);

    // image utils functions
    template<typename pixel_t>
    friend image_base<pixel_t> pad(const image_base<pixel_t> & input, std::vector<int> pre, std::vector<int> post);
    
    template<typename pixel_t>
    friend image_base<pixel_t> unpad(const image_base<pixel_t> & input, std::vector<int> pre, std::vector<int> post);

    template<typename pixel_t_in, typename pixel_t_out>
    friend image_base<pixel_t_out> cast(const image_base<pixel_t_in> & input, pixel_t_out cast_type);
    
    template<typename pixel_t>
    friend image_base<pixel_t> normalize(const image_base<pixel_t> & input, pixel_t min, pixel_t max);
    
    // friend functions
    // template<typename pixel_t>
    // friend image_base<std::complex<pixel_t>> real2complex(const image_base<pixel_t> & input);
    template<typename pixel_t>
    friend image_base<pixel_t> complex2real(const image_base<std::complex<pixel_t>> & input);
    
    // fft in 2d only now
    template<typename pixel_t>
    friend image_base<std::complex<pixel_t>> fft(const image_base<pixel_t> & input);
    template<typename pixel_t>
    friend image_base<std::complex<pixel_t>> ifft(const image_base<std::complex<pixel_t>> & input);
    
    template<typename pixel_t>
    friend typename image_base<pixel_t>::vector gradient(const image_base<pixel_t> & input);
    // std::vector<image_base<pixel_type>> gradient();

    // Access
    // virtual ptr_pixels2 neighbors(int e);
    virtual ptr_pixels4 neighbors4(int e);
    virtual ptr_pixels8 neighbors8(int e);

    // TODO
    // create operator << to print info of image as image_info function [DONE] (inherited object)
    // class operations: +,-,*,/  [DONE]
    // initialize data with zeros, ones, random [DONE]
    // filters: normalize (0 to 1), padding, gaussian, convolution [not], gradient, fft [DONE]
    // scalar operations: scalar*Image [DONE]
    // functions in_place: transpose, add, substract, multiply, divide, pow
    // extra functions: copy, cast [DONE]
};


// ===========================================
//      Functions of Class image_base
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base()
{
    init(0, 0);
    data.reset();
};

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(int d)
{
    assert(d>1);
    assert(d<4);
    if (d == 2){ init(0, 0); };
    if (d == 3){ init(0, 0, 0); };
    data.reset();
};

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(int w, int h)
{
    init(w, h);
    allocate(width*height);
    // data = std::shared_ptr<pixel_type[]>(new pixel_type[width*height]);
    // data = std::make_shared<pixel_type[]>(width*height);
    // data = std::make_shared<float[]>(width*height); // dynamic allocation
};

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(int w, int h, int l)
{
    init(w, h, l);
    allocate(width*height*length);
};

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(std::vector<int> size)
{
    if (size.size() == 2)
    {
        init(size[0], size[1]);
        allocate(size[0]*size[1]);
    }
    else if (size.size() == 3)
    {
        init(size[0], size[1], size[2]);
        allocate(size[0]*size[1]*size[2]);
    }
    else
    {
        ; // undefined for other dimensions
    };
};

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(ptr_vector buffer, int w, int h)
{
    init(w, h);
    data.reset();
    data = buffer;
};

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(ptr_vector buffer, int w, int h, int l)
{
    init(w, h, l);
    data.reset();
    data = buffer;
};

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(std::initializer_list<pixel_type> list)
{
    init(list.size(), 1);
    allocate(list.size());

    pixel_type * p1 = this->ptr();
    for (auto i = list.begin(); i != list.end(); i++)
    {
        *p1 = *i;
        p1++;
        // std::cout << *i << std::endl;
    };
}

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(const image_base<pixel_type> & input)
{
    copy(input);
};

// Destructor
// template <typename pixel_type>
// image_base<pixel_type>::~image_base()
// {
//     data.reset();
// };

template <typename pixel_type>
template <typename ... ARGS>
typename image_base<pixel_type>::pointer image_base<pixel_type>::new_pointer(const ARGS&... args)
{
    return std::make_shared<image_base<pixel_type>>(args...); // not working for inherited classes
};

// Empty init
template <typename pixel_type>
void image_base<pixel_type>::init(int w, int h)
{   
    this->class_name = "image_base";
    channels = 1;
    width = w;
    height = h;
    length = 1;
    num_elements = width*height;

    space_object<pixel_type>::init(2);

    this->size = std::vector<int>{width, height};
};

template <typename pixel_type>
void image_base<pixel_type>::init(int w, int h, int l)
{   
    this->class_name = "image_base";
    channels = 1;
    width = w;
    height = h;
    length = l;
    num_elements = width*height*length;

    space_object<pixel_type>::init(3);

    this->size = std::vector<int>{width, height, length};
};

// Copy metadata
template <typename pixel_type>
void image_base<pixel_type>::copy_properties(const image_base<pixel_type> & input)
{
    width = input.get_width();
    height = input.get_height();
    length = input.get_length();
    channels = input.get_channels();
    num_elements = input.get_total_elements();

    space_object<pixel_type>::mimic_(input);
};

// Allocate
template <typename pixel_type>
void image_base<pixel_type>::allocate(int total_elements)
{
    data.reset();
    data = std::make_shared<std::vector<pixel_type>>(total_elements);
};

// template <typename pixel_type>
// typename image_base<pixel_type>::pointer image_base<pixel_type>::anew()
// {
//     return std::make_shared<image_base<pixel_type>>();
// };

// Full copy
template <typename pixel_type>
void image_base<pixel_type>::clone(const image_base<pixel_type> & input)
{
    copy_properties(input);
    allocate(input.get_total_elements());

    pixel_type * p1 = this->ptr();
    pixel_type * p2 = input.ptr();

    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p1[k] = p2[k];
    };
};

// Point to the same data
template <typename pixel_type>
void image_base<pixel_type>::copy(const image_base<pixel_type> & input)
{
    copy_properties(input);
    data.reset();
    data = input.get_data();
};

template <typename pixel_type>
void image_base<pixel_type>::mimic(const image_base<pixel_type> & input)
{
    copy_properties(input);
    data.reset();
    allocate(input.get_total_elements());
};

// template <typename pixel_type>
// typename image_base<pixel_type>::pointer image_base<pixel_type>::new_pointer()
// {
//     image_base<pixel_type>::pointer img = std::make_shared<image_base<pixel_type>>();
//     return img;
// };


// ===========================================
// Interface Functions
// ===========================================
template <typename pixel_type>
void image_base<pixel_type>::read(std::string file_name)
{
    if (this->get_dimension() == 2) { read_2d(file_name); };
    if (this->get_dimension() == 3) { read_3d(file_name); };
};

template <typename pixel_type>
void image_base<pixel_type>::read_2d(std::string file_name)
{
    const int d = 2;
    using itkImageType = itk::Image<pixel_type, d>;
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
    pixel_type * p = image_itk->GetBufferPointer();
    int w = itksize[0];
    int h = itksize[1];

    init(w, h);
    data.reset();
    allocate(num_elements);

    // Copy of the image
    // #pragma omp parallel for // ** Better times without omp
    for(int k=0; k<num_elements; k++)
    {
        (*data)[k] = *(p+k);
    };

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

template <typename pixel_type>
void image_base<pixel_type>::read_3d(std::string file_name)
{
    const int d = 3;
    using itkImageType = itk::Image<pixel_type, d>;
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
    pixel_type * p = image_itk->GetBufferPointer();
    int w = itksize[0];
    int h = itksize[1];
    int l = itksize[2];

    init(w, h, l);
    data.reset();
    allocate(num_elements);

    // Copy of the image
    // #pragma omp parallel for // ** Better times without omp
    for(int k=0; k<num_elements; k++)
    {
        (*data)[k] = *(p+k);
    };

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
// void image_base::read_itk(itk::Image< type_itk, 2 > image_itk)
// {
//     template <size_t type_itk> 
//     using ImageType = itk::Image< type_itk, 2 >;
//     // typedef itk::Image< typename type_itk, 2 > ImageType;

//     ImageType::RegionType region = image_itk->GetLargestPossibleRegion();
//     ImageType::SizeType size = region.GetSize();
//     type * p = image_itk->GetBufferPointer();
//     int w = size[0];
//     int h = size[0];

//     image_base(w,h); // empty current image and create new

//     // Copy of the image
//     for(int k=0; k<w*h; k++)
//     {
//         *(data+k) = *(p+k);
//     };
// };

template <typename pixel_type>
void image_base<pixel_type>::write(std::string file_name)
{
    if (this->get_dimension() == 2) { write_2d(file_name); };
    if (this->get_dimension() == 3) { write_3d(file_name); };
};

template <typename pixel_type>
void image_base<pixel_type>::write_2d(std::string file_name)
{
    const int d = 2;
    using itkImageType = itk::Image<pixel_type, d>;
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

    pixel_type * p = image_itk->GetBufferPointer();

    for(int k=0; k<num_elements; k++)
    {
        *(p+k) = (*data)[k];
    };
    
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

template <typename pixel_type>
void image_base<pixel_type>::write_3d(std::string file_name)
{
    const int d = 3;
    using itkImageType = itk::Image<pixel_type, d>;
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

    pixel_type * p = image_itk->GetBufferPointer();

    for(int k=0; k<num_elements; k++)
    {
        *(p+k) = (*data)[k];
    };
    
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

// ===========================================
// Get Functions
// ===========================================
template <typename pixel_type>
int image_base<pixel_type>::get_width() const
{
    return width;
};

template <typename pixel_type>
int image_base<pixel_type>::get_height() const
{
    return height;
};

template <typename pixel_type>
int image_base<pixel_type>::get_length() const
{
    return length;
};

template <typename pixel_type>
int image_base<pixel_type>::get_channels() const
{
    return channels;
};

template <typename pixel_type>
int image_base<pixel_type>::get_total_elements() const
{
    return num_elements;
};

template <typename pixel_type>
int image_base<pixel_type>::get_ptr_count() const
{
    return data.use_count();
};

template <typename pixel_type>
typename image_base<pixel_type>::ptr_vector image_base<pixel_type>::get_data() const
{
    return data;
};

template <typename pixel_type>
pixel_type * image_base<pixel_type>::ptr() const
{
    if (num_elements > 0)  { return data.get()->data(); }
    else                   { return nullptr; };
    // return (*data).data();
};

template <typename pixel_type>
typename image_base<pixel_type>::iterator image_base<pixel_type>::begin() const
{
    return data.get()->begin();
};

template <typename pixel_type>
typename image_base<pixel_type>::iterator image_base<pixel_type>::end() const
{
    return data.get()->end();
};

// template <typename pixel_type>
// std::shared_ptr<std::vector<image_base<pixel_type>>> image_base<pixel_type>::get_grid()
// {
//     return grid;
// };

// ===========================================
// Print Functions
// ===========================================
template <typename pixel_type>
void image_base<pixel_type>::print_ptr_count()
{
    std::cout << "internal image ptr count: ";
    std::cout << data.use_count() << std::endl; // print existing shared pointer
};

template <typename pixel_type>
std::string image_base<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Image Information";
    if (msg != "") { title = msg; };
    // Summary of the image information
    ss << space_object<pixel_type>::info(title);
    // ss << "\n===== " << title << " =====\n";
    // ss << "Pixel type: \t\t" << typeid((*data)[0]).name() << std::endl;
    // ss << "Dimensions: \t\t" << dim << std::endl;
    ss << "Pixel channels: \t" << channels << std::endl; 
    
    // ss << "Size: \t\t\t[ ";
    // for(int i = 0; i < this->size.size(); i++) { ss << this->size[i] << " "; };
    // ss << "]" << std::endl;
    
    // ss << "Length (mm): \t\t[ ";
    // for(int i = 0; i < this->spacing.size(); i++) { ss << this->spacing[i] << " "; };
    // ss << "]" << std::endl;

    // ss << "Origin (mm): \t\t[ ";
    // for(int i = 0; i < this->origin.size(); i++) { ss << this->origin[i] << " "; };
    // ss << "]" << std::endl;

    ss << "Total elements: \t" << num_elements << std::endl; //Get the total number of pixels

    ss << "Data pointer: \t\t" << ptr() << std::endl; //Get the total number of pixels
    // ss << std::endl;

    return ss.str();
};

template <typename pixel_type>
std::string image_base<pixel_type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; };

    if (this->dim == 2)
    {
        pixel_type * p = this->ptr();
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
        int w = this->width;
        int h = this->height;
        int l = this->length;
        pixel_type * p = this->ptr();

        // std::cout << "Image data:" << std::endl;
        // ss << "[";
        for(int k = 0; k < l; k++)
        {
            ss << "[ ";
            for(int i = 0; i < h; i++)
            {
                for(int j=0; j < w; j++)
                {
                    ss << p[j + i*w + k*w*h] << " "; // valgrind error solved
                };
                if(i < h-1){ss << std::endl << "  ";};
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
template <typename pixel_type>
void image_base<pixel_type>::zeros()
{
    pixel_type * p = this->ptr();
    // std::cout << num_elements << std::endl;
    
    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p[k] = (pixel_type)0; // casting to pixel_type
    };
};

template <typename pixel_type>
void image_base<pixel_type>::ones()
{
    pixel_type * p = this->ptr();

    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p[k] = (pixel_type)1.0; // casting to pixel_type
    };
};

template <typename pixel_type>
void image_base<pixel_type>::random(float min, float max)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    // std::default_random_engine gen(rd()); //Standard random generator()
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> uniform(min, max);
    pixel_type * p = this->ptr();

    // if (typeid(pixel_type) == typeid(std::complex))
    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p[k] = (pixel_type)uniform(gen); // casting to pixel_type
    };
        
};

// Set all pixel to a fixed value
template <typename pixel_type>
void image_base<pixel_type>::fill(pixel_type value)
{
    pixel_type * p = this->ptr();

    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p[k] = value; // casting to pixel_type
    };
};

// ===========================================
// Overloading Operators
// ===========================================
// Access
template <typename pixel_type>
pixel_type & image_base<pixel_type>::operator () (int e)
{
    // assert(e < num_elements);
    return (*data)[e];
};

template <typename pixel_type>
pixel_type & image_base<pixel_type>::operator () (int w, int h)
{
    // assert(w < width and h < height);
    return (*data)[w+h*width];
};

template <typename pixel_type>
pixel_type & image_base<pixel_type>::operator () (int w, int h, int l)
{
    // assert(w < width and h < height and l < length);
    return (*data)[w+h*width+l*height*width];
};

template <typename pixel_type>
void image_base<pixel_type>::assert_size(const image_base<pixel_type> & input)
{
    assert(this->get_width()  == input.get_width());
    assert(this->get_height() == input.get_height());
    assert(this->get_length() == input.get_length());
    assert(this->get_channels() == input.get_channels());
    return;
};

// Equal
template <typename pixel_type>
image_base<pixel_type> & image_base<pixel_type>::operator = (const image_base<pixel_type> & input)
{
    // delete &data;
    copy(input);
    return *this;
};

// Image to Image
template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator + (const image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert_size(input);

    // int w = input.get_width();
    // int h = input.get_height();
    // int elements = w*h;
    int elements = input.get_total_elements();
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(input);

    // Create pointers
    // std::cout << "create correct\n";
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = input.ptr();
    pixel_type * p3 = result.ptr();

    // std::cout << "pointers correct\n";

    // #pragma omp parallel for
    for(int k=0; k<elements; k++)
    {
        p3[k] = p1[k] + p2[k];
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator - (const image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert_size(input);

    // int w = input.get_width();
    // int h = input.get_height();
    // int elements = w*h;
    int elements = input.get_total_elements();
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(input);

    // Create pointers
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = input.ptr();
    pixel_type * p3 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<elements; k++)
    {
        p3[k] = p1[k] - p2[k];
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator * (const image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert_size(input);

    // int w = input.get_width();
    // int h = input.get_height();
    // int elements = w*h;
    int elements = input.get_total_elements();
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(input);

    // Create pointers
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = input.ptr();
    pixel_type * p3 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<elements; k++)
    {
        p3[k] = p1[k] * p2[k];
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator / (const image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert_size(input);

    // int w = input.get_width();
    // int h = input.get_height();
    // int elements = w*h;
    int elements = input.get_total_elements();
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(input);

    // Create pointers
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = input.ptr();
    pixel_type * p3 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<elements; k++)
    {
        p3[k] = p1[k] / p2[k];
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator ^ (const image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert_size(input);

    // int w = input.get_width();
    // int h = input.get_height();
    // int elements = w*h;
    int elements = input.get_total_elements();
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(input);

    // Create pointers
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = input.ptr();
    pixel_type * p3 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<elements; k++)
    {
        p3[k] = pow( p1[k] , p2[k] );
    };
    return result;
};

// Scalar right hand side
template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator + (pixel_type scalar)
{
    image_base<pixel_type> result(*this); // init a image with same poperties
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p2[k] = p1[k] + scalar;
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator - (pixel_type scalar)
{
    image_base<pixel_type> result(*this);
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p2[k] = p1[k] - scalar;
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator * (pixel_type scalar)
{
    image_base<pixel_type> result(*this);
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p2[k] = p1[k] * scalar;
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator / (pixel_type scalar)
{
    image_base<pixel_type> result(*this);
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p2[k] = p1[k] / scalar;
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator ^ (pixel_type scalar)
{
    image_base<pixel_type> result(*this);
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p2[k] = pow( p1[k] , scalar );
    };
    return result;
};

// Scalar left hand side defined in header due to use of friend functions
template <typename pixel_type>
image_base<pixel_type> operator + (pixel_type scalar, image_base<pixel_type> & input)
{
    return input + scalar;
};

template <typename pixel_type>
image_base<pixel_type> operator - (pixel_type scalar, const image_base<pixel_type> & input)
{
    image_base<pixel_type> result(input);
    pixel_type * p1 = input.ptr();
    pixel_type * p2 = result.ptr();
    // std::cout << "elements: " << input.get_total_elements() << "\n";

    // #pragma omp parallel for
    for(int k=0; k<input.get_total_elements(); k++)
    {
        p2[k] = scalar - p1[k];
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> operator * (pixel_type scalar, image_base<pixel_type> & input)
{
    return input * scalar;
};

template <typename pixel_type>
image_base<pixel_type> operator / (pixel_type scalar, const image_base<pixel_type> & input)
{
    image_base<pixel_type> result(input);
    pixel_type * p1 = input.ptr();
    pixel_type * p2 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<input.get_total_elements(); k++)
    {
        p2[k] = scalar / p1[k];
    };
    return result;
};

// ===========================================
// Functions
// ===========================================
// Matrix product
template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::_x_(const image_base<pixel_type> & input)
{
    using Map = Eigen::Map<Eigen::Matrix<pixel_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
    image_base<pixel_type> result(input.width, height);
    // std::shared_ptr<image_base<pixel_type>> result = 
    // std::shared_ptr<image_base<pixel_type>>(new image_base<pixel_type>(height, input.width));
    // result.print();
    // std::cout << "ptr: " << result.get_data() << std::endl;

    Map A(ptr(), height, width);
    Map B(input.ptr(), input.get_height(), input.get_width());
    Map C(result.ptr(), result.get_height(), result.get_width());

    // std::cout << A << std::endl;
    // std::cout << std::endl;
    // std::cout << B << std::endl;
    // std::cout << std::endl;
    // std::cout << A*B << std::endl;
    // std::cout << std::endl;

    C = (A*B);
    // C = (A.transpose()*B.transpose()).transpose(); // transpose function doesn't cost anything (average=12ns)
    return result;
};

// Reductions operations
// Reductions not using OpenMP. They are slower.
template <typename pixel_type>
pixel_type image_base<pixel_type>::min()
{
    pixel_type x = 0;
    pixel_type * p1 = this->ptr();
    if (num_elements > 0) { x = p1[0]; };

    // #pragma omp parallel for reduction(min:x)
    for(int k=1; k<num_elements; k++)
    {
        x = std::min(x,p1[k]);
    };
    return x;
};

template <typename pixel_type>
pixel_type image_base<pixel_type>::max()
{
    pixel_type x = 0;
    pixel_type * p1 = this->ptr();
    if (num_elements > 0) { x = p1[0]; };

    // #pragma omp parallel for reduction(max:x)
    for(int k=1; k<num_elements; k++)
    {
        x = std::max(x,p1[k]);
    };
    return x;
};

template <typename pixel_type>
pixel_type image_base<pixel_type>::sum()
{
    pixel_type x = 0;
    pixel_type * p1 = this->ptr();
    // std::cout << "elements: " << num_elements << std::endl;

    // #pragma omp parallel for reduction(+:x)
    for(int k=0; k<num_elements; k++)
    {
        x += p1[k];
    };
    return x;
};

template <typename pixel_type>
pixel_type image_base<pixel_type>::prod()
{
    pixel_type x = 1;
    pixel_type * p1 = this->ptr();

    // #pragma omp parallel for reduction(*:x)
    for(int k=0; k<num_elements; k++)
    {
        x *= p1[k];
    };
    return x;
};

// Vectorial dot product. Verify the same number of elements, then product and reduce
template <typename pixel_type>
pixel_type image_base<pixel_type>::dot(const image_base<pixel_type> & input)
{
    // check same number of elements
    assert(num_elements == input.get_total_elements());

    // using the builtin functions
    // pixel_type x;
    // image_base<pixel_type> result = (*this)*input;
    // x = result.sum();

    // from scratch
    pixel_type x = 0;
    pixel_type * p1 = this->ptr();
    pixel_type * p2 = input.ptr();

    for(int k=0; k<num_elements; k++)
    {
        x += p1[k]*p2[k];
    };

    return x;
};

template<typename pixel_type>
image_base<pixel_type> pad(const image_base<pixel_type> & input, std::vector<int> pre, std::vector<int> post)
{
    int w = input.get_width();
    int h = input.get_height();
    int l = input.get_length();
    int N = input.get_total_elements();
    pixel_type * p1 = input.ptr();

    std::vector<int> extra(pre.size(),0);
    for (int i = 0; i < extra.size(); i++){ extra[i] = pre[i]+post[i]; };
    
    image_base<pixel_type> result;
    if (input.get_dimension() == 2){ result = image_base<pixel_type>(w+extra[0], h+extra[1]); };
    if (input.get_dimension() == 3){ result = image_base<pixel_type>(w+extra[0], h+extra[1], l+extra[2]); };
    result.zeros();
    pixel_type * p2 = result.ptr();

    int ww = result.get_width();
    int hh = result.get_height();
    int ll = result.get_length();
    int c = 0;
    if (input.get_dimension() == 2)
    {
        for(int i = pre[1]; i < hh-post[1]; i++)
        {
            for(int j = pre[0]; j < ww-post[0]; j++)
            {
                p2[j + i*ww] = p1[c];
                c++;
            };
        };
    };
    if (input.get_dimension() == 3)
    {
        for (int k = pre[2]; k < ll-post[2]; k++)
        {
            for(int i = pre[1]; i < hh-post[1]; i++)
            {
                for(int j = pre[0]; j < ww-post[0]; j++)
                {
                    p2[j + i*ww + k*ww*hh] = p1[c];
                    c++;
                };
            };
        };
    };
    
    return result;
};

template<typename pixel_type>
image_base<pixel_type> unpad(const image_base<pixel_type> & input, std::vector<int> pre, std::vector<int> post)
{
    int w = input.get_width();
    int h = input.get_height();
    int l = input.get_length();
    int N = input.get_total_elements();
    pixel_type * p1 = input.ptr();

    std::vector<int> extra(pre.size(),0);
    for (int i = 0; i < extra.size(); i++){ extra[i] = pre[i]+post[i]; };
    
    image_base<pixel_type> result;
    if (input.get_dimension() == 2){ result = image_base<pixel_type>(w-extra[0], h-extra[1]); };
    if (input.get_dimension() == 3){ result = image_base<pixel_type>(w-extra[0], h-extra[1], l-extra[2]); };
    result.zeros();
    pixel_type * p2 = result.ptr();

    int c = 0;
    if (input.get_dimension() == 2)
    {
        for(int i = pre[1]; i < h-post[1]; i++)
        {
            for(int j = pre[0]; j < w-post[0]; j++)
            {
                p2[c] = p1[j + i*w];
                c++;
            };
        };
    };
    if (input.get_dimension() == 3)
    {
        for (int k = pre[2]; k < l-post[2]; k++)
        {
            for(int i = pre[1]; i < h-post[1]; i++)
            {
                for(int j = pre[0]; j < w-post[0]; j++)
                {
                    p2[c] = p1[j + i*w + k*w*h];
                    c++;
                };
            };
        };
    };
    
    return result;
};

template<typename pixel_t_in, typename pixel_t_out>
image_base<pixel_t_out> cast(const image_base<pixel_t_in> & input, pixel_t_out cast_type)
{
    image_base<pixel_t_out> result(input.get_size());
    // image_base<pixel_t_out> result;
    // if (input.get_dimension() == 2){ result = image_base<pixel_t_out>(input.get_width(), input.get_height()); };
    // if (input.get_dimension() == 3){ result = image_base<pixel_t_out>(input.get_width(), input.get_height(), input.get_length()); };
    
    pixel_t_in * p1 = input.ptr();
    pixel_t_out * p2 = result.ptr();

    // #pragma omp parallel for
    for(int k=0; k<input.get_total_elements(); k++)
    {
        p2[k] = (pixel_t_out)p1[k];
    };
    return result;
};

template<typename pixel_type>
image_base<pixel_type> normalize(const image_base<pixel_type> & input, pixel_type min = 0.0, pixel_type max = 1.0)
{
    // pixel_type minv = input.min();
    // pixel_type maxv = input.max();

    image_base<pixel_type> in;
    in.copy(input); 

    pixel_type minv = in.min();
    pixel_type maxv = in.max();

    image_base<pixel_type> result = (in - minv)*((max - min)/(maxv - minv)) + min;
    return result;
};

template<typename pixel_type>
image_base<pixel_type> complex2real(const image_base<std::complex<pixel_type>> & input)
{
    int d = input.get_dimension();
    int N = input.get_total_elements();
    std::complex<pixel_type> * p1 = input.ptr();

    image_base<pixel_type> result(input.get_size());
    
    pixel_type * p2 = result.ptr();

    for (int i = 0; i < N; i++)
    {
        p2[i] = std::real(p1[i]);
    };
    return result;
};

template <typename pixel_type>
image_base<std::complex<pixel_type>> fft(const image_base<pixel_type> & input)
{
    int w = input.get_width();
    int h = input.get_height();
    int N = input.get_total_elements();
    fftw_complex in[N], out[N];
    pixel_type * p1 = input.ptr();

    for (int i = 0; i < N; i++)
    {
        in[i][0] = p1[i];
        in[i][1] = 0;
    };

    fftw_plan p_fft;
    p_fft = fftw_plan_dft_2d(h, w, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p_fft);
    fftw_destroy_plan(p_fft);

    image_base<std::complex<pixel_type>> result(w,h);
    std::complex<pixel_type> * p2 = result.ptr();    
    for (int i = 0; i < N; i++)
    {
        p2[i] = std::complex(out[i][0],out[i][1]);
    };

    return result;
};

template <typename pixel_type>
image_base<std::complex<pixel_type>> ifft(const image_base<std::complex<pixel_type>> & input)
{
    int w = input.get_width();
    int h = input.get_height();
    int N = input.get_total_elements();
    fftw_complex in[N], out[N];
    std::complex<pixel_type> * p1 = input.ptr();

    std::complex<pixel_type> v;
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

    image_base<std::complex<pixel_type>> result(w,h);
    std::complex<pixel_type> * p2 = result.ptr();    
    for (int i = 0; i < N; i++)
    {
        p2[i] = std::complex(out[i][0]/N,out[i][1]/N);
    };

    return result;
};

template<typename pixel_type>
typename image_base<pixel_type>::vector gradient(const image_base<pixel_type> & input)
{
    std::vector<int> none{0,0};    
    std::vector<int> extra{1,1};
    image_base<pixel_type> input_pad = pad(input, none, extra);
    // input_pad.print_data();
    int w = input_pad.get_width();
    int h = input_pad.get_height();

    // derivative x
    image_base<pixel_type> dx(w,h);
    pixel_type * pdx = dx.ptr();
    pdx[0] = 0.5;
    pdx[1] = 0;
    pdx[2] = -0.5;
    // dx.print_data();

    // derivative y
    image_base<pixel_type> dy(w,h);
    pixel_type * pdy = dy.ptr();
    pdy[0] = 0.5;
    pdy[w] = 0.0;
    pdy[2*w] = -0.5;
    // dy.print_data();

    typename image_base<pixel_type>::vector grad(2);
    typename image_base<pixel_type>::pointer didx(new image_base<pixel_type>(w,h));
    typename image_base<pixel_type>::pointer didy(new image_base<pixel_type>(w,h));
    grad[0] = didx;
    grad[1] = didy;

    // image_base<pixel_type> didx(w,h);
    image_base<pixel_type> didx_ = complex2real(ifft(fft(input_pad)*fft(dx)));
    *grad[0] = unpad(didx_, std::vector<int>{1,0}, std::vector<int>{0,1});
    // image_base<pixel_type> didx = unpad(didx_, std::vector<int>{1,0}, std::vector<int>{0,1});

    // image_base<pixel_type> didy(w,h);
    image_base<pixel_type> didy_ = complex2real(ifft(fft(input_pad)*fft(dy)));
    *grad[1] = unpad(didy_, std::vector<int>{0,1}, std::vector<int>{1,0});
    // image_base<pixel_type> didy = unpad(didy_, std::vector<int>{0,1}, std::vector<int>{1,0});

    // return didx;
    // return didy;
    return grad;
};

template <typename pixel_type>
typename image_base<pixel_type>::ptr_pixels4 image_base<pixel_type>::neighbors4(int e)
{
    ptr_pixels4 arr = std::make_unique<std::array<pixel_type,4>>();
    // std::cout << "image_base";
    return arr;
}

template <typename pixel_type>
typename image_base<pixel_type>::ptr_pixels8 image_base<pixel_type>::neighbors8(int e)
{
    ptr_pixels8 arr = std::make_unique<std::array<pixel_type,8>>();
    return arr;
}

}; //end namespace

// Template constructions
// template class image_base<unsigned char>;  // 1 byte
// template class image_base<unsigned short>; // 2 byte
// template class image_base<short>;          // 2 byte
// template class image_base<unsigned int>;   // 4 byte
// template class image_base<int>;            // 4 byte
// template class image_base<float>;          // 4 byte
// template class image_base<double>;         // 8 byte

#endif