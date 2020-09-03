/*
* @Author: jose
* @Date:   2020-08-24 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-08-24 00:00:00
*/

#ifndef __CUDA_BUFFER_H__
#define __CUDA_BUFFER_H__

// std libs
#include <iostream>     // std::cout
#include <vector>       // std::vector

// local libs
#include "object.h"
#include "cuda/interface.cuh"
#include "cuda/kernels.cuh"

namespace imart
{

template <typename type>
class cuda_buffer: public inherit<cuda_buffer<type>, object>
{
public:
    // Type definitions
    using self    = cuda_buffer;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    type * x;
    int _size_;

    // ===========================================
    // Functions
    // ===========================================
    void init(int sz);
    void create_memory();
    void clean_memory();

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    cuda_buffer();
    // Constructors
    cuda_buffer(int sz);
    // Destructor
    ~cuda_buffer();

    // ===========================================
    // Get Functions
    // ===========================================
    type * get() const;
    int size() const;

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    void push_memory(type * data, int sz = -1, int offset = 0);
    void push_memory(const type * data, int sz = -1, int offset = 0);
    void pull_memory(type * data, int sz = -1, int offset = 0);
};

// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
// Constructor
template <typename type>
cuda_buffer<type>::cuda_buffer()
{   
    this->class_name = "cuda_buffer";
    init(0);
};

// Constructor
template <typename type>
cuda_buffer<type>::cuda_buffer(int sz)
{
    this->class_name = "cuda_buffer";
    init(sz);
};

template <typename type>
cuda_buffer<type>::~cuda_buffer()
{
    clean_memory();
};

// ===========================================
// Create Functions
// ===========================================
template <typename type>
void cuda_buffer<type>::init(int sz)
{
    // Set internal variables
    _size_ = sz;
    // Allocate
    create_memory();
    // if (_size_ > 0)
    //     create_memory();
    // else x = nullptr;
    // std::cout << "create buffer" << std::endl;
};

template <typename type>
void cuda_buffer<type>::create_memory()
{
    cuda_create_memory(x, _size_);
};

template <typename type>
void cuda_buffer<type>::clean_memory()
{
    cuda_clean_memory(x);
};

// ===========================================
// Get Functions
// ===========================================
template <typename type>
type * cuda_buffer<type>::get() const
{
    // if (_size_ == 0) return nullptr;
    return x;
};

template <typename type>
int cuda_buffer<type>::size() const
{
    return _size_;
};

// ===========================================
// Set Functions
// ===========================================

// ===========================================
// Print Functions
// ===========================================
template <typename type>
std::string cuda_buffer<type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "CUDA buffer Information";
    if (msg != "") { title = msg; };
    // Summary of the optimizer information
    ss << object::info(title);

    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type>
void cuda_buffer<type>::push_memory(type * data, int sz, int offset)
{
    if (sz == -1)
        cuda_push_memory(x, data, _size_, offset);
    else if (sz > 0)
        cuda_push_memory(x, data, sz, offset);
    else ;
};

template <typename type>
void cuda_buffer<type>::push_memory(const type * data, int sz, int offset)
{
    if (sz == -1)
        cuda_push_memory(x, data, _size_, offset);
    else if (sz > 0)
        cuda_push_memory(x, data, sz, offset);
    else ;
};

template <typename type>
void cuda_buffer<type>::pull_memory(type * data, int sz, int offset)
{
    if (sz == -1)
        cuda_pull_memory(x, data, _size_, offset);
    else if (sz > 0)
        cuda_pull_memory(x, data, sz, offset);
    else ;
};

}; //end namespace

#endif