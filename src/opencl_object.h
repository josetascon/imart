/*
* @Author: jose
* @Date:   2020-01-27 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-01-27 00:00:00
*/

#ifndef __OPENCL_OBJECT_H__
#define __OPENCL_OBJECT_H__


// std libs
#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <cassert>      // assert

// gpu libs
#include <CL/cl.hpp>
#include "object.h"

namespace imart
{

class opencl_object: public inherit<opencl_object, object>
{
public:
    //Type definitions
    ;
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    cl_int _err_;
    cl::Device _device_;
    cl::Platform _platform_;
    cl::Context _context_;
    cl::Program _program_;
    cl::CommandQueue _queue_;
    cl::Kernel _kernel_;
    cl_int work_group_size;
    int count;

    bool status_init;
    bool status_program;
    bool status_kernel;

    // ===========================================
    // Functions
    // ===========================================
    void init();
    void check_error(int err);
    template<typename Arg, typename ...Args>
    void argument(Arg const & first, Args const &... args);
    // template<typename Arg>
    // void argument(Arg&& arg);
    // template<typename ...Args>
    // void argument(Args&&... args);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    opencl_object();

    // ===========================================
    // Get Functions
    // ===========================================
    cl::Platform get_platform() const;
    cl::Device get_device() const;
    cl::Context get_context() const;
    cl::Program get_program() const;
    cl::CommandQueue get_queue() const;
    cl::Kernel get_kernel() const;
    cl_int get_work_group_size() const;

    // ===========================================
    // Set Functions
    // ===========================================
    // int set_device(cl::Device device);
    // int set_context(cl::Context context);
    void set_program(std::string code);
    void set_kernel(std::string kernel_name);

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // compute
    void program(std::string code, std::string kernel_name);
    template<typename ...Args>
    void arguments(Args&&... args);
    void execute(int max_size);
};

// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
// Constructor
opencl_object::opencl_object()
{   
    this->class_name = "opencl_object";
    init();
};

void opencl_object::init()
{
    // Set internal variables
    status_init = false;
    status_program = false;
    status_kernel = false;
    _err_ = 0;

    // Get the platform
    std::vector<cl::Platform> platforms;
    _err_ = cl::Platform::get(&platforms);
    assert(platforms.size() > 0);
    _platform_ = platforms.front();

    // Get the devices in the computer
    std::vector<cl::Device> devices;
    _err_ = _platform_.getDevices(CL_DEVICE_TYPE_GPU, &devices);
    assert(devices.size() > 0);
    _device_ = devices.front();

    // Create the context and the queue
    _context_ = cl::Context(_device_);
    _queue_ = cl::CommandQueue(_context_, _device_);

    // Initialization succeded
    status_init = true;
};

// ===========================================
// Get Functions
// ===========================================
cl::Platform opencl_object::get_platform() const
{
    assert(status_init);
    return _platform_;
};

cl::Device opencl_object::get_device() const
{
    assert(status_init);
    return _device_;
};

cl::Context opencl_object::get_context() const
{
    assert(status_init);
    return _context_;
};

cl::CommandQueue opencl_object::get_queue() const
{
    assert(status_init);
    return _queue_;
};

cl::Program opencl_object::get_program() const
{
    assert(status_program);
    return _program_;
};

cl::Kernel opencl_object::get_kernel() const
{
    assert(status_kernel);
    return _kernel_;
};

cl_int opencl_object::get_work_group_size() const
{
    return work_group_size;
};

// ===========================================
// Set Functions
// ===========================================
void opencl_object::set_program(std::string code)
{
    // Reset error value
    _err_ = 0;

    // Read the kernel in a external file    
    cl::Program::Sources sources(1, std::make_pair(code.c_str(), code.length() + 1));

    // Create the program
    _program_ = cl::Program(_context_, sources);
    _err_ = _program_.build("-cl-std=CL1.2");
    assert(_err_ == 0);

    // Created program
    status_program = true;
};

void opencl_object::set_kernel(std::string kernel_name)
{
    // Reset error value
    _err_ = 0;

    // Add kernel
    _kernel_ = cl::Kernel(_program_, kernel_name.c_str(), &_err_);
    assert(_err_ == 0);

    // Read work group size
    work_group_size = _kernel_.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(_device_, &_err_);
    // work_group_size = 32;
    assert(_err_ == 0);

    // Created program
    status_kernel = true;
};

void opencl_object::program(std::string code, std::string kernel_name)
{
    set_program(code);
    set_kernel(kernel_name);
};

template<typename Arg, typename ...Args>
void opencl_object::argument(Arg const & first, Args const &... args)
{
    // std::cout << "Set Kernel Argument: " << count << std::endl;
    _kernel_.setArg(count, first);
    count = count + 1;
    if constexpr (sizeof...(args) > 0) 
    {
        argument(args...);
    };
};

template<typename ...Args>
void opencl_object::arguments(Args&&... args)
{
    count = 0;
    argument(args...);
};

void opencl_object::execute(int max_size)
{
    if (max_size < work_group_size) work_group_size = max_size;
    // std::cout << "Run Kernel" << std::endl;
    // _err_ = _queue_.enqueueNDRangeKernel(_kernel_, cl::NullRange, cl::NDRange(max_size), cl::NDRange(work_group_size));
    _err_ = _queue_.enqueueNDRangeKernel(_kernel_, cl::NullRange, cl::NDRange(max_size));
    // std::cout << "Kernel error" << _err_ << std::endl;
    assert(_err_ == 0);
};

std::string opencl_object::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "OpenCL Object Information";
    if (msg != "") { title = msg; };
    // Summary of the optimizer information
    ss << object::info(title);
    // if(status_init)
    // if(status_program)
    // if(status_kernel)
    return ss.str();
};

}; //end namespace

#endif