/*
* @Author: jose
* @Date:   2019-01-21 09:36:17
* @Last Modified by:   jose
* @Last Modified time: 2019-01-21 09:40:26
*/

#ifndef __OPENCL_CONFIG_H__
#define __OPENCL_CONFIG_H__

#include <fstream>
#include <string>
#include <cassert>

#include <CL/cl.hpp>

cl::Program create_program(const std::string& file);








cl::Program create_program(const std::string& file)
{

    // Get the platform
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    assert(platforms.size() > 0);
    // std::cout << "Platform\tOK\n";

    // Get the devices in the computer
    auto platform = platforms.front();
    std::vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);

    assert(devices.size() > 0);

    auto device = devices.front();
    
    // Read the kernel in a external file
    std::ifstream hello_world_file(file);
    std::string src(std::istreambuf_iterator<char>(hello_world_file), (std::istreambuf_iterator<char>()));

    cl::Program::Sources sources(1, std::make_pair(src.c_str(), src.length() + 1));

    // Create the context and the program
    cl::Context context(device);
    cl::Program program(context, sources);

    auto err = program.build("-cl-std=CL2.1");
    return program;
}


#endif