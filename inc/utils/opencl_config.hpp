/*
* @Author: jose
* @Date:   2019-01-21 09:36:17
* @Last Modified by:   jose
* @Last Modified time: 2019-01-21 09:40:26
*/

#ifndef __OPENCL_CONFIG_HPP__
#define __OPENCL_CONFIG_HPP__

#include <fstream>
#include <string>
#include <cassert>

#include <CL/cl.hpp>

cl::Program create_program(const std::string& file);

#endif