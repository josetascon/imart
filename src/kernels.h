/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __KERNELS_H__
#define __KERNELS_H__

// std libs
#include <iostream>     // std::cout
#include <string>       // std::string
#include <typeinfo>     // operator typeids

namespace imart
{

template<typename type>
std::string string_type(){ return typeid(type).name();};
template<>
std::string string_type<unsigned int>(){ return "unsigned int";};
template<>
std::string string_type<int>(){ return "int";};
template<>
std::string string_type<float>(){ return "float";};
template<>
std::string string_type<double>(){ return "double";};

std::string kernel_scalar(std::string input_type, std::string op, bool function=false, bool reverse=false)
{
    std::string input1, input2;
    if (reverse){ input1 = "scalar"; input2 = "vin[gid]"; }
    else { input1 = "vin[gid]"; input2 = "scalar"; }

    std::string source;
    source.append("__kernel void kernel_scalar(");
    source.append("     __global const " + input_type + " * vin, ");
    source.append("     __global " + input_type + " * vout, ");
    source.append(      input_type + " scalar)\n");
    source.append("{\n");
    source.append("    int gid  = get_global_id(0);\n");    
    if(function){ source.append("    vout[gid] = " + op + "(" + input1 + " , " + input2 + ");\n"); }
    else { source.append("    vout[gid] = " + input1 + " " + op + " " + input2 + ";\n"); };
    source.append("};\n;");
    return source;
};

std::string kernel_vector(std::string input_type, std::string op, bool function=false, bool reverse=false)
{
    std::string input1, input2;
    if (reverse){ input1 = "vin2[gid]"; input2 = "vin1[gid]"; }
    else { input1 = "vin1[gid]"; input2 = "vin2[gid]"; }

    std::string source;
    source.append("__kernel void kernel_vector(");
    source.append("     __global const " + input_type + " * vin1, ");
    source.append("     __global const " + input_type + " * vin2, ");
    source.append("     __global " + input_type + " * vout)\n");
    source.append("{\n");
    source.append("    int gid  = get_global_id(0);\n");    
    if(function){ source.append("    vout[gid] = " + op + "(" + input1 + " , " + input2 + ");\n"); }
    else { source.append("    vout[gid] = " + input1 + " " + op + " " + input2 + ";\n"); };
    source.append("};\n;");
    return source;
};

std::string kernel_cast(std::string input_type, std::string output_type)
{
    std::string source;
    source.append("__kernel void kernel_cast(");
    source.append("     __global const " + input_type + " * vin, ");
    source.append("     __global " + output_type + " * vout)\n");
    source.append("{\n");
    source.append("    int gid  = get_global_id(0);\n");    
    source.append("    vout[gid] = (" + output_type + ")vin[gid];\n");
    source.append("};\n;");
    return source;
};

}; //end namespace

#endif