/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-15 13:55:00
*/


#ifndef __TRANSFORM_BASE_H__
#define __TRANSFORM_BASE_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <cassert>      // assert

// images 
#include "image_base.h"
#include "grid.h"


// extra matrix eigen
#include <eigen3/Eigen/Core>

// parallel
// openmp
// opencl


// Definitions
// typedef float pixel_type;

// Class to allow friend methods
// template <typename pixel_type>
// class transform_base;

// template <typename pixel_type>
// std::ostream & operator << (std::ostream & os, transform_base<pixel_type> & input);


// Class image_base_2d
template <typename pixel_type>
class transform_base
{
protected:
    int dim;
    std::string type;
    image_base<pixel_type> parameters;

public:
    // ===========================================
    // Create Functions
    // ===========================================
    transform_base();
    transform_base(image_base<pixel_type> & params);

    // ===========================================
    // Get Functions
    // ===========================================
    int get_dimension();
    std::string get_type();
    image_base<pixel_type> get_parameters();


    // ===========================================
    // Print Functions
    // ===========================================
    void print(std::string msg = "");
    std::string info(std::string msg = "");

    template <typename pixel_type_>
    friend std::ostream & operator << (std::ostream & os, transform_base<pixel_type_> & input);
    
    // ===========================================
    // Initialization Functions
    // ===========================================
    virtual void identity(); // parameters that do nothing when transform is applied

    // ===========================================
    // Functions
    // ===========================================
    virtual std::vector<pixel_type> transform(std::vector<pixel_type> & point);
    virtual grid<pixel_type> transform(grid<pixel_type> & input);

};


// ===========================================
//          Functions of Class transform_base
// ===========================================

// ===========================================
// Create Functions
// ===========================================
template <typename pixel_type>
transform_base<pixel_type>::transform_base()
{
    dim = 1;
    type = "transform base";
    parameters = image_base<pixel_type>();
};

template <typename pixel_type>
transform_base<pixel_type>::transform_base(image_base<pixel_type> & params)
{
    dim = 1;
    type = "transform base";
    parameters = params;
};

// ===========================================
// Print Functions
// ===========================================
template <typename pixel_type>
void transform_base<pixel_type>::print(std::string msg)
{
    std::cout << transform_base::info(msg);
};

template <typename pixel_type>
std::string transform_base<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Transform Information";
    if (msg != "") { title = msg; };

    // Summary of the image information
    ss << "\n===== " << title << " =====\n";
    
    ss << "Pixel type: \t\t" << typeid(parameters(0)).name() << std::endl;
    ss << "Transform type: \t" << type << std::endl; 
    ss << "Dimensions: \t\t" << dim << std::endl;
    
    ss << "Parameters: \t\t";
    ss << parameters.info_data();
    // ss << "]" << std::endl;
    ss << std::endl;

    return ss.str();
};


template <typename pixel_type>
std::ostream & operator << (std::ostream & os, transform_base<pixel_type> & input)
{
    os << input.info("");
    return os;
};


// ===========================================
// Initialization Functions
// ===========================================
template <typename pixel_type>
void transform_base<pixel_type>::identity()
{
    ;
};

// ===========================================
// Functions
// ===========================================
template <typename pixel_type>
std::vector<pixel_type> transform_base<pixel_type>::transform(std::vector<pixel_type> & point)
{
    return point;
};

template <typename pixel_type>
grid<pixel_type> transform_base<pixel_type>::transform(grid<pixel_type> & input)
{
    return input;
};

#endif