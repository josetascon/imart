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
#include <typeinfo>     // operator typeid
#include <cassert>      // assert
#include <cmath>        // math (pow)

// images itk
#include <itkImage.h>
#include <itkImageFileReader.h>

// extra matrix eigen
#include <eigen3/Eigen/Core>

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
class image_base
{
    //Type definitions
    // template<typename pixel_t> using vector = std::vector<image_base<pixel_t>>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int dim;
    int width;
    int height;
    int length;
    int num_elements;
    int channels;       // for now will only support one channel

    std::vector<int> size;
    std::vector<pixel_type> spacing;
    std::vector<pixel_type> origin;
    std::vector<pixel_type> direction;

    // Image data
    std::shared_ptr<pixel_type[]> data; // change float to template

    // std::shared_ptr<vector_image> grid;

public:

    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    //! Constructor empty.
    image_base();
    //! Constructor using width and height.
    image_base(int w, int h); //ONLY 2d***
    //! Constructor with existing data.
    image_base(std::shared_ptr<pixel_type[]> buffer, int w, int h); //ONLY 2d***
    //! Constructor with file_name (call read()).
    image_base(std::string file_name);
    //! Destructor
    ~image_base();

    // consider these methods as protected
    void init(int w, int h); //ONLY 2d***
    void update(image_base<pixel_type> & input);

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
    //! Get the image dimension.
    int get_dimension();
    //! Get the image width.
    int get_width();
    //! Get the image height.
    int get_height();
    //! Get the image length.
    int get_length();
    //! Get the number of elements allocated.
    int get_total_elements();
    
    std::vector<int> get_size();
    std::vector<pixel_type> get_spacing();
    std::vector<pixel_type> get_origin();
    std::vector<pixel_type> get_direction();
    std::shared_ptr<pixel_type[]> get_data();

    int get_ptr_count();

    // std::shared_ptr<vector_image> get_grid();

    // ===========================================
    // Print Functions
    // ===========================================
    void print(std::string msg = "");
    void print_data(std::string msg = "");
    void print_ptr_count();
    std::string info(std::string msg = "");
    std::string info_data(std::string msg = ""); //ONLY 2d***

    template<typename pixel_t>
    friend std::ostream & operator << (std::ostream & os, image_base<pixel_t> & input);

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
    pixel_type operator () (int e);
    pixel_type operator () (int w, int h); //ONLY 2d***

    // Equal
    image_base<pixel_type> & operator = (image_base<pixel_type> input);
    
    // Image to Image
    image_base<pixel_type> operator + (image_base<pixel_type> & input);
    image_base<pixel_type> operator - (image_base<pixel_type> & input);
    image_base<pixel_type> operator * (image_base<pixel_type> & input);
    image_base<pixel_type> operator / (image_base<pixel_type> & input);
    image_base<pixel_type> operator ^ (image_base<pixel_type> & input);


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
    friend image_base<pixel_t> operator - (pixel_t scalar, image_base<pixel_t> & input);
    template<typename pixel_t>
    friend image_base<pixel_t> operator * (pixel_t scalar, image_base<pixel_t> & input);
    template<typename pixel_t>
    friend image_base<pixel_t> operator / (pixel_t scalar, image_base<pixel_t> & input);

    // ===========================================
    // Functions
    // ===========================================
    // Matrix product //ONLY 2d***
    image_base<pixel_type> _x_(image_base<pixel_type> & input);
    // virtual grid meshgrid(); // moved to grid class

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
    data = std::shared_ptr<pixel_type[]>(new pixel_type[width*height]);
    // data = std::make_shared<float[]>(width*height); // dynamic allocation
};

// Constructor
template <typename pixel_type>
image_base<pixel_type>::image_base(std::shared_ptr<pixel_type[]> buffer, int w, int h)
{
    init(w, h);
    data.reset();
    data = buffer;
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
    dim = 2;
    channels = 1;
    width = w;
    height = h;
    num_elements = width*height;

    size = std::vector<int>{width, height};
    spacing = std::vector<pixel_type>(dim, 1.0);
    origin = std::vector<pixel_type>(dim, 0.0);
    direction = std::vector<pixel_type>(dim*dim);

    // initialize direction, identity matrix
    int den = dim + 1;
    for(int i=0; i < dim*dim; i++){ if((i%den)==0) { direction[i] = 1.0; }; };
};

// Copy metadata
template <typename pixel_type>
void image_base<pixel_type>::update(image_base<pixel_type> & input)
{
    width = input.get_width();
    height = input.get_height();
    length = input.get_length();
    num_elements = input.get_total_elements();

    size = input.get_size();
    spacing = input.get_spacing();
    origin = input.get_origin();
    direction = input.get_direction();
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
    typename itkImageType::SizeType size = region.GetSize();
    pixel_type * p = image_itk->GetBufferPointer();
    int w = size[0];
    int h = size[1];

    // image_base(w,h); // empty current image and create new

    this->init(w, h);
    data.reset();
    data = std::shared_ptr<pixel_type[]>(new pixel_type[width*height]);

    // Copy of the image
    #pragma omp parallel for
    for(int k=0; k<w*h; k++)
    {
        data[k] = *(p+k);
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
int image_base<pixel_type>::get_dimension()
{
    return dim;
};


template <typename pixel_type>
int image_base<pixel_type>::get_width()
{
    return width;
};

template <typename pixel_type>
int image_base<pixel_type>::get_height()
{
    return height;
};

template <typename pixel_type>
int image_base<pixel_type>::get_length()
{
    return length;
};

template <typename pixel_type>
int image_base<pixel_type>::get_total_elements()
{
    return num_elements;
};

template <typename pixel_type>
std::vector<int> image_base<pixel_type>::get_size()
{
    return size;
};

template <typename pixel_type>
std::vector<pixel_type> image_base<pixel_type>::get_spacing()
{
    return spacing;
};

template <typename pixel_type>
std::vector<pixel_type> image_base<pixel_type>::get_origin()
{
    return origin;
};

template <typename pixel_type>
std::vector<pixel_type> image_base<pixel_type>::get_direction()
{
    return direction;
};

template <typename pixel_type>
std::shared_ptr<pixel_type[]> image_base<pixel_type>::get_data()
{
    return data;
};

template <typename pixel_type>
int image_base<pixel_type>::get_ptr_count()
{
    return data.use_count();
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
void image_base<pixel_type>::print(std::string msg)
{
    std::cout << this->info(msg);
};

template <typename pixel_type>
void image_base<pixel_type>::print_data(std::string msg)
{
    std::cout << this->info_data(msg) << std::endl;
};

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
    ss << "\n===== " << title << " =====\n";
    
    ss << "Pixel type: \t\t" << typeid(data[0]).name() << std::endl;
    ss << "Pixel channels: \t" << channels << std::endl; 
    ss << "Dimensions: \t\t" << dim << std::endl;
    
    ss << "Size: \t\t\t[ ";
    for(int i = 0; i < size.size(); i++) { ss << size[i] << " "; };
    ss << "]" << std::endl;
    
    ss << "Length (mm): \t\t[ ";
    for(int i = 0; i < spacing.size(); i++) { ss << spacing[i] << " "; };
    ss << "]" << std::endl;

    ss << "Origin (mm): \t\t[ ";
    for(int i = 0; i < origin.size(); i++) { ss << origin[i] << " "; };
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
    // std::cout << "Image data:" << std::endl;
    // std::cout << "["
    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            ss << data[j+i*width] << " ";
        };
        ss << std::endl;
    };
    // std::cout << "]";
    // ss << std::endl;
    return ss.str();
};

template <typename pixel_type>
std::ostream & operator << (std::ostream & os, image_base<pixel_type> & input)
{
    os << input.info("");
    return os;
};

// ===========================================
// Initialization Functions
// ===========================================
template <typename pixel_type>
void image_base<pixel_type>::zeros()
{
    #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        data[k] = (pixel_type)0.0; // casting to pixel_type
    };
};

template <typename pixel_type>
void image_base<pixel_type>::ones()
{
    #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        data[k] = (pixel_type)1.0; // casting to pixel_type
    };
};

template <typename pixel_type>
void image_base<pixel_type>::random(pixel_type min, pixel_type max)
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    // std::default_random_engine gen(rd()); //Standard random generator()
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> uniform(min, max);

    #pragma omp parallel for
    for(int k=0; k<num_elements; k++)
    {
        data[k] = (pixel_type)uniform(gen); // casting to pixel_type
    };
        
};

// ===========================================
// Overloading Operators
// ===========================================
// Access
template <typename pixel_type>
pixel_type image_base<pixel_type>::operator () (int e)
{
    // assert(e < num_elements);
    return this->get_data()[e];
};

template <typename pixel_type>
pixel_type image_base<pixel_type>::operator () (int w, int h)
{
    // assert(w < width && h < height);
    return this->get_data()[w+h*width];
};

// Equal
template <typename pixel_type>
image_base<pixel_type> & image_base<pixel_type>::operator = (image_base<pixel_type> input)
{
    // delete &data;
    this->data.reset();
    this->data = input.get_data();
    this->update(input);
    return *this;
};

// Image to Image
template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator + (image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

    int w = input.get_width();
    int h = input.get_height();
    int elements = w*h;
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(w, h);

    // Create pointers
    // std::cout << "create correct\n";
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = input.get_data();
    std::shared_ptr<pixel_type[]> p3 = result.get_data();

    // std::cout << "pointers correct\n";

    #pragma omp parallel for
    for(int k=0; k<elements; k++)
    {
        // p1[k] = p1[k] + p2[k];
        p3[k] = p1[k] + p2[k];
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator - (image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

    int w = input.get_width();
    int h = input.get_height();
    int elements = w*h;
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(w, h);

    // Create pointers
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = input.get_data();
    std::shared_ptr<pixel_type[]> p3 = result.get_data();

    #pragma omp parallel for
    for(int k=0; k<elements; k++)
    {
        p3[k] = p1[k] - p2[k];
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator * (image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

    int w = input.get_width();
    int h = input.get_height();
    int elements = w*h;
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(w, h);

    // Create pointers
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = input.get_data();
    std::shared_ptr<pixel_type[]> p3 = result.get_data();

    #pragma omp parallel for
    for(int k=0; k<elements; k++)
    {
        p3[k] = p1[k] * p2[k];
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator / (image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

    int w = input.get_width();
    int h = input.get_height();
    int elements = w*h;
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(w, h);

    // Create pointers
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = input.get_data();
    std::shared_ptr<pixel_type[]> p3 = result.get_data();

    #pragma omp parallel for
    for(int k=0; k<elements; k++)
    {
        p3[k] = p1[k] / p2[k];
    };
    return result;
};

template <typename pixel_type>
image_base<pixel_type> image_base<pixel_type>::operator ^ (image_base<pixel_type> & input)
{
    // std::cout << "assert images size" << std::endl;
    assert(this->get_width()==input.get_width());
    assert(this->get_height()==input.get_height());

    int w = input.get_width();
    int h = input.get_height();
    int elements = w*h;
    // std::cout << "Image size: [" << w << ", "<< h << "]\n";
    image_base<pixel_type> result(w, h);

    // Create pointers
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = input.get_data();
    std::shared_ptr<pixel_type[]> p3 = result.get_data();

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
    image_base<pixel_type> result(width, height);
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = result.get_data();

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
    image_base<pixel_type> result(width, height);
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = result.get_data();

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
    image_base<pixel_type> result(width, height);
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = result.get_data();

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
    image_base<pixel_type> result(width, height);
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = result.get_data();

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
    image_base<pixel_type> result(width, height);
    std::shared_ptr<pixel_type[]> p1 = this->get_data();
    std::shared_ptr<pixel_type[]> p2 = result.get_data();

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
image_base<pixel_type> operator - (pixel_type scalar, image_base<pixel_type> & input)
{
    static image_base<pixel_type> result(input.get_width(), input.get_height());
    std::shared_ptr<pixel_type[]> p1 = input.get_data();
    std::shared_ptr<pixel_type[]> p2 = result.get_data();

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
image_base<pixel_type> operator / (pixel_type scalar, image_base<pixel_type> & input)
{
    static image_base<pixel_type> result(input.get_width(), input.get_height());
    std::shared_ptr<pixel_type[]> p1 = input.get_data();
    std::shared_ptr<pixel_type[]> p2 = result.get_data();

    #pragma omp parallel for
    for(int k=0; k<input.get_total_elements(); k++)
    {
        p2[k] = scalar/p1[k];
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

    Map A(get_data().get(), height, width);
    Map B(input.get_data().get(), input.get_height(), input.get_width());
    Map C(result.get_data().get(), result.get_height(), result.get_width());

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

// moved to grid class
// template <typename pixel_type>
// void image_base<pixel_type>::meshgrid()
// {
//     ;
// };

// Template constructions
// template class image_base<unsigned char>;  // 1 byte
// template class image_base<unsigned short>; // 2 byte
// template class image_base<short>;          // 2 byte
// template class image_base<unsigned int>;   // 4 byte
// template class image_base<int>;            // 4 byte
// template class image_base<float>;          // 4 byte
// template class image_base<double>;         // 8 byte

#endif