/*
* @Author: jose
* @Date:   2019-11-05 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __IMAGE_BASE_HPP__
#define __IMAGE_BASE_HPP__

// std libs
#include <iostream>     // std::cout
#include <sstream>      // stringstream
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <random>       // random
#include <typeinfo>     // operator typeid
#include <assert.h>       // assert

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

// Class image_base_2d
template <typename pixel_type>
class image_base_2d
{
private:
    // ===========================================
    // Internal Variables
    // ===========================================
    int dim;
    int width;
    int height;
    int num_elements;

    int channels;       // for now will only support one channel

    std::vector<int> size;
    std::vector<pixel_type> spacing;
    std::vector<pixel_type> origin;

    // Image data
    std::shared_ptr<pixel_type[]> data; // change float to template

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    //! Constructor empty.
    image_base_2d();

    //! Constructor using width and height.
    /*!
      Constructor that allocate memory type <float> 
      using width and height
      param w an integer with width
      param h an integer with width
    */
    image_base_2d(int w, int h);

    //! Constructor with existing data.
    /*!
      Constructor that allocate memory of <type> 
      using width and height
      param buffer a shared pointer <type> where the data is allocated
      param w an integer with width
      param h an integer with width
    */
    image_base_2d(std::shared_ptr<pixel_type[]> buffer, int w, int h);

    //Destructor
    ~image_base_2d();    

    void init(int w, int h);

    void update(image_base_2d<pixel_type> & input);

    // ===========================================
    // Interface Functions
    // ===========================================
    // interface with ITK, eigen?
    void read(std::string file_name);

    void write();

    // template <size_t type_itk>
    // void read_itk(itk::Image< type_itk, 2 > image_itk);
    // void write_itk();

    // ===========================================
    // Get Functions
    // ===========================================
    //! Get the image width.
    /*!
      \return int value with width.
    */
    int get_width();

    //! Get the image height.
    /*!
      \return int value with height.
    */
    int get_height();

    //! Get the number of elements allocated.
    /*!
      \return int value with the number of elements (pixels)
    */
    int get_total_elements();
    
    std::vector<int> get_size();

    std::vector<pixel_type> get_spacing();

    std::vector<pixel_type> get_origin();

    std::shared_ptr<pixel_type[]> get_data();

    int get_ptr_count();

    // ===========================================
    // Print Functions
    // ===========================================
    void print(std::string msg = "");
    void print_data(std::string msg = "");
    void print_ptr_count();
    std::string info(std::string msg = "");

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
    pixel_type operator () (int w, int h);

    // Equal
    image_base_2d<pixel_type> & operator = (image_base_2d<pixel_type> input);
    
    // Image to Image
    image_base_2d<pixel_type> operator + (image_base_2d<pixel_type> & input);
    image_base_2d<pixel_type> operator - (image_base_2d<pixel_type> & input);
    image_base_2d<pixel_type> operator * (image_base_2d<pixel_type> & input);
    image_base_2d<pixel_type> operator / (image_base_2d<pixel_type> & input);

    // Scalar
    image_base_2d<pixel_type> operator + (pixel_type scalar);
    image_base_2d<pixel_type> operator - (pixel_type scalar);
    image_base_2d<pixel_type> operator * (pixel_type scalar);
    image_base_2d<pixel_type> operator / (pixel_type scalar);

    // Friend classes to support double side
    friend image_base_2d<pixel_type> operator + (pixel_type scalar, image_base_2d<pixel_type> & input)
    {
        return input + scalar;
    };

    // friend image_base_2d<pixel_type> operator - (pixel_type scalar, image_base_2d<pixel_type> & input)
    // {
    //     static image_base_2d<pixel_type> result(input.get_width(), input.get_height());
    //     std::shared_ptr<pixel_type[]> p1 = input.get_data();
    //     std::shared_ptr<pixel_type[]> p2 = result.get_data();

    //     for(int k=0; k<input.num_elements; k++)
    //     {
    //         p2[k] = scalar - p1[k];
    //     };
    //     return result;
    // };

    friend image_base_2d<pixel_type> operator * (pixel_type scalar, image_base_2d<pixel_type> & input)
    {
        return input * scalar;
    };

    friend image_base_2d<pixel_type> operator / (pixel_type scalar, image_base_2d<pixel_type> & input)
    {
        static image_base_2d<pixel_type> result(input.get_width(), input.get_height());
        std::shared_ptr<pixel_type[]> p1 = input.get_data();
        std::shared_ptr<pixel_type[]> p2 = result.get_data();

        for(int k=0; k<input.num_elements; k++)
        {
            p2[k] = scalar/p1[k];
        };
        return result;
    };

    // ===========================================
    // Functions
    // ===========================================
    // Matrix product
    image_base_2d<pixel_type> _x_(image_base_2d<pixel_type> & input);

    // TODO
    // create operator << to print info of image as image_info function
    // class operations: +,-,*,/  [DONE]
    // initialize data with zeros, ones, random [DONE]
    // filters: normalize (0 to 1), padding, gaussian, convolution, gradient, fft?
    // scalar operations: scalar*Image
    // functions in_place: transpose, add, substract, multiply, divide, pow
    // extra functions: copy, cast

};

template <typename pixel_type>
inline image_base_2d<pixel_type> operator - (pixel_type scalar, image_base_2d<pixel_type> & input)
{
    static image_base_2d<pixel_type> result(input.get_width(), input.get_height());
    std::shared_ptr<pixel_type[]> p1 = input.get_data();
    std::shared_ptr<pixel_type[]> p2 = result.get_data();

    for(int k=0; k<input.get_total_elements(); k++)
    {
        p2[k] = scalar - p1[k];
    };
    return result;
};

// Template constructions
template class image_base_2d<unsigned char>;  // 1 byte
template class image_base_2d<unsigned short>; // 2 byte
template class image_base_2d<short>;          // 2 byte
template class image_base_2d<unsigned int>;   // 4 byte
template class image_base_2d<int>;            // 4 byte
template class image_base_2d<float>;          // 4 byte
template class image_base_2d<double>;         // 8 byte

#endif