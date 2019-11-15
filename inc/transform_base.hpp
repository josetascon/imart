


#ifndef __TRANSFORM_BASE_HPP__
#define __TRANSFORM_BASE_HPP__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <assert.h>       // assert

// images 
#include "image_base.hpp"

// extra matrix eigen
#include <eigen3/Eigen/Core>

// parallel
// openmp
// opencl


// Definitions
// typedef float pixel_type;

// Class image_base_2d
template <typename pixel_type>
class transform_base
{
private:
    int dim;
    std::string type;
    image_base<pixel_type> parameters;

public:
    transform_base();

    // ===========================================
    // Print Functions
    // ===========================================
    int get_dimension();
    std::string get_type();
    image_base<pixel_type> get_parameters();


    // ===========================================
    // Print Functions
    // ===========================================
    void print(std::string msg = "");
    std::string info(std::string msg = "");

    // operator << 

    // ===========================================
    // Functions
    // ===========================================
    void transform(std::vector<pixel_type> & point);
    void transform(image_base<pixel_type> & image);

};

#endif