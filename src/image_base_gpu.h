/*
* @Author: jose
* @Date:   2019-12-03 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-12-03 00:00:00
*/

#ifndef __IMAGE_BASE_GPU_H__
#define __IMAGE_BASE_GPU_H__

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
#include "image_base_gpu.h"


// Class image_base_gpu
template <typename pixel_type>
class image_base_gpu: public object<pixel_type>
{
public:
    //Type definitions
    using iterator = typename std::vector<pixel_type>::iterator;
    using ptr_vector = std::shared_ptr<std::vector<pixel_type>>;
    using ptr_pixels4 = std::unique_ptr<std::array<pixel_type,4>>;
    using ptr_pixels8 = std::unique_ptr<std::array<pixel_type,8>>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int width;
    int height;
    int length;
    int num_elements;
    int channels;       // for now will only support one channel

    // Image data
    viennacl::vector<pixel_type> data;
    // MISSIN CODE. TAKE IT FROM image_base

public:

};

#endif