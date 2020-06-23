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


std::string kernel_assign(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_assign(");
    source.append("     __global " + input_type + " * vin, ");
    source.append("     " + input_type + " value)\n");
    source.append("{\n");
    source.append("    int gid  = get_global_id(0);\n");    
    source.append("    vin[gid] = value;\n");
    source.append("};");
    return source;
};

std::string kernel_copy(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_copy(");
    source.append("     __global const " + input_type + " * vin, ");
    source.append("     __global " + input_type + " * vout)\n");
    source.append("{\n");
    source.append("    int gid  = get_global_id(0);\n");    
    source.append("    vout[gid] = vin[gid];\n");
    source.append("};");
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
    source.append("};");
    return source;
};

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
    source.append("};");
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
    source.append("};");
    return source;
};

std::string kernel_sum(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_sum(");
    source.append("     __global const " + input_type + " * vin, ");
    source.append("     __local " + input_type + " * local_data, ");
    source.append("     __global " + input_type + " * vout)\n");
    source.append("{\n");
    source.append("     size_t gid = get_global_id(0);\n");
    source.append("     size_t local_size = get_local_size(0);\n");
    source.append("     size_t lid = get_local_id(0);\n");
    source.append("     local_data[lid] = vin[gid];\n");
    source.append("     barrier(CLK_LOCAL_MEM_FENCE);\n");
    source.append("     for (unsigned int i = local_size >> 1; i > 0; i >>= 1)\n");
    source.append("     {\n");
    source.append("         if (lid < i)\n");
    source.append("         {\n");
    source.append("             local_data[lid] += local_data[lid + i];\n");
    source.append("         };\n");
    source.append("         barrier(CLK_LOCAL_MEM_FENCE);\n");
    source.append("     };\n");
    source.append("     barrier(CLK_LOCAL_MEM_FENCE);\n");
    source.append("     if(lid == 0)\n");
    source.append("     {\n");
    source.append("         vout[get_group_id(0)] = local_data[0];\n");
    source.append("     };\n");
    source.append("};");
    return source;
};

std::string kernel_minmax(std::string input_type, bool is_max)
{
    std::string source;
    if (is_max)
        source.append("__kernel void kernel_max(");
    else
        source.append("__kernel void kernel_min(");
    source.append("     __global const " + input_type + " * vin, ");
    source.append("     int len, ");
    source.append("     __local " + input_type + " * local_data, ");
    source.append("     __global " + input_type + " * vout)\n");
    source.append("{\n");
    source.append("     size_t gid = get_global_id(0);\n");
    source.append("     size_t local_size = get_local_size(0);\n");
    source.append("     size_t lid = get_local_id(0);\n");
    // source.append("     local_data[lid] = vin[gid];\n");
    // source.append("     barrier(CLK_LOCAL_MEM_FENCE);\n");
    source.append("     " + input_type + " thread_result = vin[0]; \n");
    source.append("     for (unsigned int i = get_global_id(0); i < len; i += get_global_size(0)) \n");
    source.append("     { \n");
    source.append("         " + input_type + " tmp = vin[i];\n");
    if (is_max)
        source.append("         thread_result = thread_result > tmp ? thread_result : tmp;\n");
    else
        source.append("         thread_result = thread_result < tmp ? thread_result : tmp;\n");
    source.append("     };\n");
    source.append("  local_data[get_local_id(0)] = thread_result; \n");

    source.append("     for (unsigned int i = local_size >> 1; i > 0; i >>= 1)\n");
    source.append("     {\n");
    source.append("         if (lid < i)\n");
    source.append("         {\n");
    if (is_max)
        source.append("             local_data[lid] =  local_data[lid] > local_data[lid + i]? local_data[lid] : local_data[lid + i];\n");
    else
        source.append("             local_data[lid] =  local_data[lid] < local_data[lid + i]? local_data[lid] : local_data[lid + i];\n");
    source.append("         };\n");
    source.append("         barrier(CLK_LOCAL_MEM_FENCE);\n");
    source.append("     };\n");
    source.append("     barrier(CLK_LOCAL_MEM_FENCE);\n");
    source.append("     if(lid == 0)\n");
    source.append("     {\n");
    source.append("         vout[get_group_id(0)] = local_data[0];\n");
    source.append("     };\n");
    source.append("};");
    return source;
};

std::string kernel_random(std::string input_type, float minv, float maxv, bool seed=true )
{
    // CHECK RAND function in python
    std::string source;
    source.append("float rand(float seed, float minv, float maxv)\n");
    source.append("{\n");
    source.append("     float a = (maxv-minv)/2147483947.0;\n");
    source.append("     float b = minv;\n");
    source.append("     float tmp = fmod(seed*16807.0,2147483947.0);\n");
    source.append("     return a*tmp + b;\n");
    source.append("};\n");
    // source.append("float mfunc(float * minv)\n");
    // source.append("{\n");
    // source.append("     return *minv;\n");
    // source.append("};\n");
    if (seed) source.append("__kernel void kernel_random(");
    else      source.append("__kernel void kernel_random_next(");
    source.append("     __global " + input_type + " * vin)\n");
    source.append("{\n");
    source.append("    int gid  = get_global_id(0);\n");
    // source.append("    float a = 18.0;\n");
    // source.append("    vin[gid] = ("+input_type+")rand(&a);\n");

    if (seed) source.append("    float seed = (float)gid;\n");
    else      source.append("    float seed = (float)vin[gid];\n");
    source.append("    vin[gid] = ("+input_type+")rand(seed,0.0,1.0);\n");
    // source.append("    vin[gid] = ("+input_type+")rand(1.0,");
    // source.append("    " + std::to_string(minv) + "," + std::to_string(maxv) + ");\n");
    source.append("};");
    return source;
};

std::string kernel_grid_2d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_grid_2d(");
    source.append("     __global " + input_type + " * x,");
    source.append("     __global " + input_type + " * y,");
    source.append("     __global double * sod,"); // consider conversion to float to support all gpu
    source.append("     int w,");
    source.append("     int h)\n");
    source.append("{\n");
    source.append("     double s0 = sod[0]; double s1 = sod[1];\n");
    source.append("     double o0 = sod[2]; double o1 = sod[3];\n");
    source.append("     double d0 = sod[4]; double d1 = sod[5];\n");
    source.append("     double d2 = sod[6]; double d3 = sod[7];\n");
    source.append("     for (unsigned int j = get_global_id(1); j < h; j += get_global_size(1)) \n");
    source.append("     {\n");
    source.append("         for (unsigned int i = get_global_id(0); i < w; i += get_global_size(0)) \n");
    source.append("         {\n");
    source.append("             x[i+j*w] = (" + input_type + ")d0*s0*i + d1*s1*j + o0;\n");
    source.append("             y[i+j*w] = (" + input_type + ")d2*s0*i + d3*s1*j + o1;\n");
    // source.append("             x[i+j*w] = (" + input_type + ")i;\n");
    // source.append("             y[i+j*w] = (" + input_type + ")j;\n");
    source.append("         };\n");
    source.append("     };\n");
    source.append("};");
    return source;

    // source.append("     int i  = get_global_id(0);\n");
    // source.append("     int j  = get_global_id(1);\n");
    // source.append("     if(i < w && j < h)\n");
    // source.append("     {\n");
    // source.append("         x[i+j*w] = (" + input_type + ")i;\n");
    // source.append("         y[i+j*w] = (" + input_type + ")j;\n");
    // source.append("     };\n");
};

std::string kernel_grid_3d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_grid_3d(");
    source.append("     __global " + input_type + " * x,");
    source.append("     __global " + input_type + " * y,");
    source.append("     __global " + input_type + " * z,");
    source.append("     __global double * sod,"); // consider conversion to float to support all gpu
    source.append("     int w,");
    source.append("     int h,");
    source.append("     int l)\n");
    source.append("{\n");
    source.append("     double s0 = sod[0]; double s1 = sod[1]; double s2 = sod[2];\n");
    source.append("     double o0 = sod[3]; double o1 = sod[4]; double o2 = sod[5];\n");
    source.append("     double d0 = sod[6]; double d1 = sod[7]; double d2 = sod[8];\n");
    source.append("     double d3 = sod[9]; double d4 = sod[10]; double d5 = sod[11];\n");
    source.append("     double d6 = sod[12]; double d7 = sod[13]; double d8 = sod[14];\n");
    source.append("     for (unsigned int k = get_global_id(2); k < l; k += get_global_size(2)) \n");
    source.append("     {\n");
    source.append("         for (unsigned int j = get_global_id(1); j < h; j += get_global_size(1)) \n");
    source.append("         {\n");
    source.append("             for (unsigned int i = get_global_id(0); i < w; i += get_global_size(0)) \n");
    source.append("             {\n");
    source.append("                 x[i + j*w + k*w*h] = (" + input_type + ")d0*s0*i + d1*s1*j + d2*s2*k + o0;;\n");
    source.append("                 y[i + j*w + k*w*h] = (" + input_type + ")d3*s0*i + d4*s1*j + d5*s2*k + o1;\n");
    source.append("                 z[i + j*w + k*w*h] = (" + input_type + ")d6*s0*i + d7*s1*j + d8*s2*k + o2;\n");
    // source.append("                 x[i + j*w + k*w*h] = (" + input_type + ")i;\n");
    // source.append("                 y[i + j*w + k*w*h] = (" + input_type + ")j;\n");
    // source.append("                 z[i + j*w + k*w*h] = (" + input_type + ")k;\n");
    source.append("             };\n");
    source.append("         };\n");
    source.append("     };\n");
    source.append("};");
    return source;
};

}; //end namespace

#endif