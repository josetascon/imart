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

// images itk
#include <itkImage.h>
#include <itkImageFileReader.h>

// extra matrix eigen
#include <eigen3/Eigen/Core>

// local libs
#include "object.h"

// parallel
// openmp
// opencl

// Definitions
// typedef float pixel_type;

// Class to allow friend methods
// template <typename pixel_type>
// class image_base;

// // Friend definitions
// template <typename pixel_type>
// std::ostream & operator << (std::ostream & os, image_base<pixel_type> & input);


// Class image_base
template <typename pixel_type>
class image_base: public object<pixel_type>
{
public:
    //Type definitions
    // template<typename pixel_t> using vector = std::vector<image_base<pixel_t>>;
    using iterator = typename std::vector<pixel_type>::iterator;
    using ptr_vector = std::shared_ptr<std::vector<pixel_type>>;
    using ptr_pixels4 = std::unique_ptr<std::array<pixel_type,4>>;
    using ptr_pixels8 = std::unique_ptr<std::array<pixel_type,8>>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    // int dim;
    int width;
    int height;
    int length;
    int num_elements;
    int channels;       // for now will only support one channel

    // std::vector<int> size;
    // std::vector<pixel_type> spacing;
    // std::vector<pixel_type> origin;
    // std::vector<pixel_type> direction; //not supported with interpolation

    // Image data
    ptr_vector data;
    // std::shared_ptr<pixel_type[]> data; // change float to template

    // std::shared_ptr<vector_image> grid;

    // consider these methods as protected
    virtual void init(int w, int h); //ONLY 2d***
    virtual void init(int w, int h, int l); //ONLY 3d***
    virtual void update(const image_base<pixel_type> & input); // copy only properties

public:

    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    //! Constructor empty.
    image_base();
    //! Constructor using width and height.
    image_base(int w, int h); //ONLY 2d***
    //! Constructor using width, height and length.
    image_base(int w, int h, int l); //ONLY 3d***
    //! Constructor with existing data.
    image_base(ptr_vector buffer, int w, int h); //ONLY 2d***
    //! Constructor with existing data.
    image_base(ptr_vector buffer, int w, int h, int l); //ONLY 3d***

    image_base(const image_base<pixel_type> & input);
    //! Constructor with file_name (call read()).
    image_base(std::string file_name);
    //! Destructor
    ~image_base();

    // ===========================================
    // Interface Functions
    // ===========================================
    // interface with ITK, eigen?
    void read(std::string file_name); //ONLY 2d***
    void write();

    // template <size_t type_itk>
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
    void random(pixel_type min=0.0, pixel_type max=1.0);

    // ===========================================
    // Overloading Operators
    // ===========================================
    // Access
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
    image_base<pixel_type> _x_(image_base<pixel_type> & input);
    // virtual grid meshgrid(); // moved to grid class

    // Access
    // virtual ptr_pixels2 neighbors(int e);
    virtual ptr_pixels4 neighbors4(int e);
    virtual ptr_pixels8 neighbors8(int e);


    // TODO
    // create operator << to print info of image as image_info function
    // class operations: +,-,*,/  [DONE]
    // initialize data with zeros, ones, random [DONE]
    // filters: normalize (0 to 1), padding, gaussian, convolution, gradient, fft?
    // scalar operations: scalar*Image
    // functions in_place: transpose, add, substract, multiply, divide, pow
    // extra functions: copy, cast

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
image_base<pixel_type>::image_base(int w, int h)
{
    init(w, h);
    data.reset();
    data = std::make_shared<std::vector<pixel_type>>(width*height);
    // data = std::shared_ptr<pixel_type[]>(new pixel_type[width*height]);
    // data = std::make_shared<pixel_type[]>(width*height);
    // data = std::make_shared<float[]>(width*height); // dynamic allocation
};

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(int w, int h, int l)
{
    init(w, h, l);
    data.reset();
    data = std::make_shared<std::vector<pixel_type>>(width*height*length);
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

template <typename pixel_type>
image_base<pixel_type>::image_base(const image_base<pixel_type> & input)
{
    update(input);
    // std::cout << "update\n"; 
};

// Destructor
template <typename pixel_type>
image_base<pixel_type>::~image_base()
{
    data.reset();
};

// Empty init
template <typename pixel_type>
void image_base<pixel_type>::init(int w, int h)
{   
    this->class_name = "image_base";
    channels = 1;
    width = w;
    height = h;
    num_elements = width*height;

    object<pixel_type>::init(2);
    // this->dim = 2;
    // this->size = std::vector<int>{width, height};
    // this->spacing = std::vector<pixel_type>(this->dim, 1.0);
    // this->origin = std::vector<pixel_type>(this->dim, 0.0);
    // this->direction = std::vector<pixel_type>(this->dim*this->dim);

    // // initialize direction, identity matrix
    // int den = this->dim + 1;
    // for(int i=0; i < this->dim*this->dim; i++){ if((i%den)==0) { this->direction[i] = 1.0; }; };

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

    object<pixel_type>::init(3);

    this->size = std::vector<int>{width, height, length};
};

// Copy metadata
template <typename pixel_type>
void image_base<pixel_type>::update(const image_base<pixel_type> & input)
{
    width = input.get_width();
    height = input.get_height();
    length = input.get_length();
    num_elements = input.get_total_elements();

    data.reset();
    data = std::make_shared<std::vector<pixel_type>>(num_elements);
    object<pixel_type>::update(input);
};

// ===========================================
// Interface Functions
// ===========================================
template <typename pixel_type>
void image_base<pixel_type>::read(std::string file_name)
{
    // Type definitions
    using itkImageType = itk::Image<pixel_type, 2>;
    using itkReaderType = itk::ImageFileReader<itkImageType>;
    // typedef itk::Image<pixel_type, 2>           itkImageType;
    // typedef itk::ImageFileReader<itkImageType>  itkReaderType;
    // typedef typename itkImageType::Pointer      itkImagePointerType;
    // typedef typename itkReaderType::Pointer     itkReaderPointerType;

    // // Objects
    // itkImagePointerType image_itk = itkImageType::New();
    // itkReaderPointerType reader = itkReaderType::New();
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

    // image_base(w,h); // empty current image and create new

    this->init(w, h);
    data.reset();
    data = std::make_shared<std::vector<pixel_type>>(width*height);
    // data = std::shared_ptr<std::vector<pixel_type>>(new pixel_type[width*height]);

    // Copy of the image
    // #pragma omp parallel for // ** Better times without omp
    for(int k=0; k<w*h; k++)
    {
        (*data)[k] = *(p+k);
    };

    //*** TODO COPY spacing, origin and direction
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
    return data.get()->data();
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
        // std::cout << "Image data:" << std::endl;
        // std::cout << "["
        for(int i = 0; i < height; i++)
        {
            for(int j=0; j < width; j++)
            {
                ss << (*data)[j+i*width] << " "; // valgrind error solved
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
    std::cout << num_elements << std::endl;
    
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
void image_base<pixel_type>::random(pixel_type min, pixel_type max)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    // std::default_random_engine gen(rd()); //Standard random generator()
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> uniform(min, max);
    pixel_type * p = this->ptr();

    // #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        p[k] = (pixel_type)uniform(gen); // casting to pixel_type
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

// Equal
template <typename pixel_type>
image_base<pixel_type> & image_base<pixel_type>::operator = (const image_base<pixel_type> & input)
{
    // delete &data;
    update(input);
    data.reset();
    data = input.get_data();
    return *this;
};

// Image to Image
template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator + (const image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

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

    #pragma omp parallel for
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
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

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

    #pragma omp parallel for
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
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

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

    #pragma omp parallel for
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
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

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

    #pragma omp parallel for
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
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

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

    #pragma omp parallel for
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

    #pragma omp parallel for
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

    #pragma omp parallel for
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

    #pragma omp parallel for
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

    #pragma omp parallel for
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

    #pragma omp parallel for
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

    #pragma omp parallel for
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

    #pragma omp parallel for
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
image_base<pixel_type> image_base<pixel_type>::_x_(image_base<pixel_type> & input)
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


// Template constructions
// template class image_base<unsigned char>;  // 1 byte
// template class image_base<unsigned short>; // 2 byte
// template class image_base<short>;          // 2 byte
// template class image_base<unsigned int>;   // 4 byte
// template class image_base<int>;            // 4 byte
// template class image_base<float>;          // 4 byte
// template class image_base<double>;         // 8 byte

#endif