/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __IMAGE_UTILS_H__
#define __IMAGE_UTILS_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <cassert>      // assert

// local libs
#include "image.h"
#include "grid.h"

namespace imart
{

// ===========================================
//      Functions of Vector Image
// ===========================================

template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator + (std::vector< std::shared_ptr< image<type,container> >> input1, std::vector< std::shared_ptr< image<type,container> >> & input2)
{
    assert(input1.size() == input2.size());
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = *(input1[i]) + *(input2[i]);
    }
    return output;
};

template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator - (std::vector< std::shared_ptr< image<type,container> >> input1, std::vector< std::shared_ptr< image<type,container> >> & input2)
{
    assert(input1.size() == input2.size());
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = *(input1[i]) - *(input2[i]);
    }
    return output;
};

template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator * (std::vector< std::shared_ptr< image<type,container> >> input1, std::vector< std::shared_ptr< image<type,container> >> & input2)
{
    assert(input1.size() == input2.size());
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = *(input1[i]) * (*(input2[i]));
    }
    return output;
};

template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator / (std::vector< std::shared_ptr< image<type,container> >> input1, std::vector< std::shared_ptr< image<type,container> >> & input2)
{
    assert(input1.size() == input2.size());
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = *(input1[i]) / (*(input2[i]));
    }
};

template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator * (std::vector< std::shared_ptr< image<type,container> >> input1, type scalar)
{
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = (*(input1[i])) * scalar;
    }
    return output;
};

// ===========================================
//      Functions of Gaussian
// ===========================================
template <typename type, typename container>
std::shared_ptr< image<type,container> > gaussian_kernel(int dim = 2, type sigma = 1.0, int width = 3)
{
    assert(width >= 3 and width%2 == 1);
    // Gaussian kernel is small, thats why is computed on cpu
    std::vector<int> sz(dim);
    std::vector<double> origin(dim);
    for(int i = 0; i < dim; i++)
    {
        sz[i] = width;
        origin[i] = -std::floor(width/2.0);
    };
    auto output = image<type,container>::new_pointer(sz);
    output->set_origin(origin);
    auto x = grid<type,container>::new_pointer(output);
    // x->print_data();

    if(dim == 2)
    {
        auto img = ((*x)[0]->operator^(2.0) + (*x)[1]->operator^(2.0))/((type)-2.0*pow(sigma,2.0));
        img = ((type)std::exp(1.0)) ^ ( img );
        *output = img/(img.sum());
    }
    else if(dim == 3)
    {
        auto img = ((*x)[0]->operator^(2.0) + (*x)[1]->operator^(2.0) + (*x)[2]->operator^(2.0))/((type)-2.0*pow(sigma,2.0));
        img = ((type)std::exp(1.0)) ^ ( img );
        *output = img/(img.sum());
    }

    return output;
};

// template<typename type, typename container>
// typename std::shared_ptr<image<type,container>> gaussian_filter(std::shared_ptr<image<type,container>> input, type sigma = 1.0, int kwidth = 3)
// {
//     int d = input->get_dimension();
//     auto output = input->mimic();

//     auto gkernel = gaussian_kernel<type,vector_cpu<type>>(d,sigma,kwidth);

//     if (input->get_data()->get_name() == "vector_cpu")
//     {
//         vector_cpu<type>::convolution(input->get_data(), gkernel->get_data(),
//                            output->get_data(), input->get_size(), kwidth);
//     }
//     else if (input->get_data()->get_name() == "vector_ocl")
//     {
//         // create kernel in gpu (copy from host ram)
//         auto kernel = vector_ocl<type>::new_pointer(gkernel->get_data()->size());
//         kernel->read_ram(gkernel->get_data()->data(),gkernel->get_data()->size());

//         vector_ocl<type>::convolution(input->get_data(), kernel,
//                                output->get_data(), input->get_size(), kwidth);
//     }
//     else ;
//     return output;
// };

template<typename type, typename container>
typename std::shared_ptr<image<type,container>> gaussian_filter(std::shared_ptr<image<type,container>> input, type sigma = 1.0, int kwidth = 3)
{
    int d = input->get_dimension();
    auto output = input->mimic();
    output->zeros();

    auto gkernel = gaussian_kernel<type,container>(d,sigma,kwidth);

    container::convolution(input->get_data(), gkernel->get_data(),
                           output->get_data(), input->get_size(), kwidth);
    return output;
};

}; //end namespace

#endif