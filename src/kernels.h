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
std::string string_type<short>(){ return "short";};
template<>
std::string string_type<unsigned short>(){ return "unsigned short";};
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
    source.append("    int i = get_global_id(0);\n");    
    source.append("    vin[i] = value;\n");
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
    source.append("    int i = get_global_id(0);\n");    
    source.append("    vout[i] = vin[i];\n");
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
    source.append("    int i = get_global_id(0);\n");    
    source.append("    vout[i] = (" + output_type + ")vin[i];\n");
    source.append("};");
    return source;
};

std::string kernel_scalar(std::string input_type, std::string op, bool function=false, bool reverse=false)
{
    std::string input1, input2;
    if (reverse){ input1 = "scalar"; input2 = "vin[i]"; }
    else { input1 = "vin[i]"; input2 = "scalar"; }

    std::string source;
    source.append("__kernel void kernel_scalar(");
    source.append("     __global const " + input_type + " * vin, ");
    source.append("     __global " + input_type + " * vout, ");
    source.append(      input_type + " scalar)\n");
    source.append("{\n");
    source.append("    int i = get_global_id(0);\n");    
    if(function){ source.append("    vout[i] = " + op + "(" + input1 + " , " + input2 + ");\n"); }
    else { source.append("    vout[i] = " + input1 + " " + op + " " + input2 + ";\n"); };
    source.append("};");
    return source;
};

std::string kernel_vector(std::string input_type, std::string op, bool function=false, bool reverse=false)
{
    std::string input1, input2;
    if (reverse){ input1 = "vin2[i]"; input2 = "vin1[i]"; }
    else { input1 = "vin1[i]"; input2 = "vin2[i]"; }

    std::string source;
    source.append("__kernel void kernel_vector(");
    source.append("     __global const " + input_type + " * vin1, ");
    source.append("     __global const " + input_type + " * vin2, ");
    source.append("     __global " + input_type + " * vout)\n");
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");    
    if(function){ source.append("    vout[i] = " + op + "(" + input1 + " , " + input2 + ");\n"); }
    else { source.append("    vout[i] = " + input1 + " " + op + " " + input2 + ";\n"); };
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
    source.append("     for (int i = local_size >> 1; i > 0; i >>= 1)\n");
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
    source.append("     for (int i = get_global_id(0); i < len; i += get_global_size(0)) \n");
    source.append("     { \n");
    source.append("         " + input_type + " tmp = vin[i];\n");
    if (is_max)
        source.append("         thread_result = thread_result > tmp ? thread_result : tmp;\n");
    else
        source.append("         thread_result = thread_result < tmp ? thread_result : tmp;\n");
    source.append("     };\n");
    source.append("  local_data[get_local_id(0)] = thread_result; \n");

    source.append("     for (int i = local_size >> 1; i > 0; i >>= 1)\n");
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

std::string kernel_pad_2d(std::string input_type, bool pad)
{
    std::string source;
    source.append("__kernel void kernel_pad_2d(");
    source.append("     __global const " + input_type + " * vin,");
    source.append("     __global " + input_type + " * vout,");
    source.append("     int start0,");
    source.append("     int start1,");
    source.append("     int end0,");
    source.append("     int end1)\n");
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int w = get_global_size(0);\n"); //pad with vin size, unpad vout size
    source.append("     int h = get_global_size(1);\n"); //pad with vin size, unpad vout size
    source.append("     int wo = w+start0+end0;\n");
    if (pad) source.append("     vout[start0+i+(start1+j)*wo] = vin[i+j*w];\n");
    else     source.append("     vout[i+j*w] = vin[start0+i+(start1+j)*wo];\n");
    source.append("};");
    return source;
};

std::string kernel_pad_3d(std::string input_type, bool pad)
{
    std::string source;
    source.append("__kernel void kernel_pad_3d(");
    source.append("     __global const " + input_type + " * vin,");
    source.append("     __global " + input_type + " * vout,");
    source.append("     int start0,");
    source.append("     int start1,");
    source.append("     int start2,");
    source.append("     int end0,");
    source.append("     int end1,");
    source.append("     int end2)\n");
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int k = get_global_id(2);\n");
    source.append("     int w = get_global_size(0);\n"); //vin size
    source.append("     int h = get_global_size(1);\n"); //vin size
    source.append("     int l = get_global_size(1);\n"); //vin size
    source.append("     int wo = w+start0+end0;\n"); //vout size
    source.append("     int ho = h+start1+end1;\n"); //vout size
    if (pad) source.append("     vout[start0+i+(start1+j)*wo+(start2+k)*wo*ho] = vin[i+j*w+k*w*h];\n");
    else     source.append("     vout[i+j*w+k*w*h] = vin[start0+i+(start1+j)*wo+(start2+k)*wo*ho];\n");
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
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     double c0 = sod[0]; double c1 = sod[1];\n");
    source.append("     double o0 = sod[2]; double o1 = sod[3];\n");
    source.append("     double d0 = sod[4]; double d1 = sod[5];\n");
    source.append("     double d2 = sod[6]; double d3 = sod[7];\n");
    source.append("     x[i+j*w] = (" + input_type + ")(d0*c0*i + d1*c1*j + o0);\n");
    source.append("     y[i+j*w] = (" + input_type + ")(d2*c0*i + d3*c1*j + o1);\n");
    // source.append("     x[i+j*w] = (" + input_type + ")i;\n");
    // source.append("     y[i+j*w] = (" + input_type + ")j;\n");
    source.append("};");
    return source;
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
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int k = get_global_id(2);\n");
    source.append("     double c0 = sod[0]; double c1 = sod[1]; double c2 = sod[2];\n");
    source.append("     double o0 = sod[3]; double o1 = sod[4]; double o2 = sod[5];\n");
    source.append("     double d0 = sod[6]; double d1 = sod[7]; double d2 = sod[8];\n");
    source.append("     double d3 = sod[9]; double d4 = sod[10]; double d5 = sod[11];\n");
    source.append("     double d6 = sod[12]; double d7 = sod[13]; double d8 = sod[14];\n");
    source.append("     x[i + j*w + k*w*h] = (" + input_type + ")(d0*c0*i + d1*c1*j + d2*c2*k + o0);\n");
    source.append("     y[i + j*w + k*w*h] = (" + input_type + ")(d3*c0*i + d4*c1*j + d5*c2*k + o1);\n");
    source.append("     z[i + j*w + k*w*h] = (" + input_type + ")(d6*c0*i + d7*c1*j + d8*c2*k + o2);\n");
    // source.append("     x[i + j*w + k*w*h] = (" + input_type + ")i;\n");
    // source.append("     y[i + j*w + k*w*h] = (" + input_type + ")j;\n");
    // source.append("     z[i + j*w + k*w*h] = (" + input_type + ")k;\n");
    source.append("};");
    return source;
};

std::string kernel_affine_2d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_affine_2d(");
    source.append("     __global const " + input_type + " * xin,");
    source.append("     __global const " + input_type + " * yin,");
    source.append("     __global " + input_type + " * xout,");
    source.append("     __global " + input_type + " * yout,");
    source.append("     __global const " + input_type + " * param)\n");
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n"); // Set NDRange one dimension, buffer in and out xy equal size
    source.append("     " + input_type + " a0 = param[0]; " + input_type + " a1 = param[1];\n");
    source.append("     " + input_type + " a2 = param[2]; " + input_type + " a3 = param[3];\n");
    source.append("     " + input_type + " t0 = param[4]; " + input_type + " t1 = param[5];\n");
    source.append("     xout[i] = (" + input_type + ")(a0*xin[i] + a1*yin[i] + t0);\n");
    source.append("     yout[i] = (" + input_type + ")(a2*xin[i] + a3*yin[i] + t1);\n");
    source.append("};");
    return source;
};

std::string kernel_affine_3d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_affine_3d(");
    source.append("     __global const " + input_type + " * xin,");
    source.append("     __global const " + input_type + " * yin,");
    source.append("     __global const " + input_type + " * zin,");
    source.append("     __global " + input_type + " * xout,");
    source.append("     __global " + input_type + " * yout,");
    source.append("     __global " + input_type + " * zout,");
    source.append("     __global const " + input_type + " * param)\n");
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n"); // Set NDRange one dimension, buffer in and out xyz equal size
    source.append("     " + input_type + " a0 = param[0]; " + input_type + " a1 = param[1]; " + input_type + " a2 = param[2];\n");
    source.append("     " + input_type + " a3 = param[3]; " + input_type + " a4 = param[4]; " + input_type + " a5 = param[5];\n");
    source.append("     " + input_type + " a6 = param[6]; " + input_type + " a7 = param[7]; " + input_type + " a8 = param[8];\n");
    source.append("     " + input_type + " t0 = param[9]; " + input_type + " t1 = param[10]; " + input_type + " t2 = param[11];\n");
    source.append("     xout[i] = (" + input_type + ")(a0*xin[i] + a1*yin[i] + a2*zin[i] + t0);\n");
    source.append("     yout[i] = (" + input_type + ")(a3*xin[i] + a4*yin[i] + a5*zin[i] + t1);\n");
    source.append("     zout[i] = (" + input_type + ")(a6*xin[i] + a7*yin[i] + a8*zin[i] + t2);\n");
    source.append("};");
    return source;
};

std::string kernel_affine_sod_2d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_affine_sod_2d(");
    source.append("     __global const " + input_type + " * xin,");
    source.append("     __global const " + input_type + " * yin,");
    source.append("     __global " + input_type + " * xout,");
    source.append("     __global " + input_type + " * yout,");
    source.append("     __global const double * sod)\n"); // consider conversion to float to support all gpu
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n"); // Set NDRange one dimension, buffer in and out xy equal size
    source.append("     double c0 = sod[0]; double c1 = sod[1];\n");
    source.append("     double o0 = sod[2]; double o1 = sod[3];\n");
    source.append("     double d0 = sod[4]; double d1 = sod[5];\n");
    source.append("     double d2 = sod[6]; double d3 = sod[7];\n");
    source.append("     xout[i] = (" + input_type + ")(d0*c0*xin[i] + d1*c1*yin[i] + o0);\n");
    source.append("     yout[i] = (" + input_type + ")(d2*c0*xin[i] + d3*c1*yin[i] + o1);\n");
    source.append("};");
    return source;
};

std::string kernel_affine_sod_3d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_affine_sod_3d(");
    source.append("     __global const " + input_type + " * xin,");
    source.append("     __global const " + input_type + " * yin,");
    source.append("     __global const " + input_type + " * zin,");
    source.append("     __global " + input_type + " * xout,");
    source.append("     __global " + input_type + " * yout,");
    source.append("     __global " + input_type + " * zout,");
    source.append("     __global const double * sod)\n"); // consider conversion to float to support all gpu
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n"); // Set NDRange one dimension, buffer in and out xyz equal size
    source.append("     double c0 = sod[0]; double c1 = sod[1]; double c2 = sod[2];\n");
    source.append("     double o0 = sod[3]; double o1 = sod[4]; double o2 = sod[5];\n");
    source.append("     double d0 = sod[6]; double d1 = sod[7]; double d2 = sod[8];\n");
    source.append("     double d3 = sod[9]; double d4 = sod[10]; double d5 = sod[11];\n");
    source.append("     double d6 = sod[12]; double d7 = sod[13]; double d8 = sod[14];\n");
    source.append("     xout[i] = (" + input_type + ")(d0*c0*xin[i] + d1*c1*yin[i] + d2*c2*zin[i] + o0);\n");
    source.append("     yout[i] = (" + input_type + ")(d3*c0*xin[i] + d4*c1*yin[i] + d5*c2*zin[i] + o1);\n");
    source.append("     zout[i] = (" + input_type + ")(d6*c0*xin[i] + d7*c1*yin[i] + d8*c2*zin[i] + o2);\n");
    source.append("};");
    return source;
};

std::string kernel_dfield_2d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_dfield_2d(");
    source.append("     __global const " + input_type + " * xin,");
    source.append("     __global const " + input_type + " * yin,");
    source.append("     __global const " + input_type + " * x,");
    source.append("     __global const " + input_type + " * y,");
    source.append("     __global " + input_type + " * xout,");
    source.append("     __global " + input_type + " * yout)\n");
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n"); // Set NDRange one dimension, buffer in and out xy equal size
    source.append("     xout[i] = xin[i] + x[i]);\n");
    source.append("     yout[i] = yin[i] + y[i]);\n");
    source.append("};");
    return source;
};

std::string kernel_dfield_3d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_dfield_2d(");
    source.append("     __global const " + input_type + " * xin,"); // grid coordinates
    source.append("     __global const " + input_type + " * yin,");
    source.append("     __global const " + input_type + " * zin,");
    source.append("     __global const " + input_type + " * x,");   // vector field
    source.append("     __global const " + input_type + " * y,");
    source.append("     __global const " + input_type + " * z,");
    source.append("     __global " + input_type + " * xout,");      // output coordinates
    source.append("     __global " + input_type + " * yout,");
    source.append("     __global " + input_type + " * zout)\n");
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n"); // Set NDRange one dimension, buffer in and out xy equal size
    source.append("     xout[i] = xin[i] + x[i]);\n");
    source.append("     yout[i] = yin[i] + y[i]);\n");
    source.append("     zout[i] = zin[i] + z[i]);\n");
    source.append("};");
    return source;
};

std::string kernel_nearest_interpolation_2d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_nearest_interpolation_2d(");
    source.append("     __global const " + input_type + " * xo,");
    source.append("     __global const " + input_type + " * yo,");
    source.append("     __global const " + input_type + " * imgr,");
    source.append("     __global " + input_type + " * imgo,");
    source.append("     int w,");   //img ref width
    source.append("     int h)\n"); //img ref height
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int n0 = get_global_size(0);\n");
    source.append("     int n1 = get_global_size(1);\n");

    source.append("     int x = round(xo[i + j*n0]);\n");
    source.append("     int y = round(yo[i + j*n0]);\n");
    source.append("     if(x >= 0 && x < w && y >= 0 && y < h)\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0] = imgr[x + y*w];\n");
    source.append("     };\n");
    source.append("};");
    return source;
};

std::string kernel_nearest_interpolation_3d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_nearest_interpolation_3d(");
    source.append("     __global const " + input_type + " * xo,");
    source.append("     __global const " + input_type + " * yo,");
    source.append("     __global const " + input_type + " * zo,");
    source.append("     __global const " + input_type + " * imgr,");
    source.append("     __global " + input_type + " * imgo,");
    source.append("     int w,");   //img ref width
    source.append("     int h,"); //img ref height
    source.append("     int l)\n"); //img ref length
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int k = get_global_id(2);\n");
    source.append("     int n0 = get_global_size(0);\n");
    source.append("     int n1 = get_global_size(1);\n");
    source.append("     int n2 = get_global_size(2);\n");

    source.append("     int x = round(xo[i + j*n0 + k*n0*n1]);\n");
    source.append("     int y = round(yo[i + j*n0 + k*n0*n1]);\n");
    source.append("     int z = round(yo[i + j*n0 + k*n0*n1]);\n");
    source.append("     if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l)\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = imgr[x + y*w + z*w*h];\n");
    source.append("     };\n");
    source.append("};");
    return source;
};

std::string kernel_linear_interpolation_2d(std::string input_type)
{
    std::string source;
    // source.append(input_type + " bilinear(" + input_type + "4 iv, " + input_type + " dx, " + input_type + " dy)\n");
    // source.append("{\n");
    // source.append("     " + input_type + " dxdy = dx*dy;\n");
    // source.append("     " + input_type + " r = iv.s0*(1-dx-dy+dxdy) + iv.s1*(dx-dxdy) + iv.s2*(dy-dxdy) + iv.s3*dxdy;\n");
    // source.append("     return r;\n");
    // source.append("};\n");

    source.append("__kernel void kernel_linear_interpolation_2d(");
    source.append("     __global const " + input_type + " * xo,");
    source.append("     __global const " + input_type + " * yo,");
    source.append("     __global const " + input_type + " * imgr,");
    source.append("     __global " + input_type + " * imgo,");
    source.append("     int w,");   //img ref width
    source.append("     int h)\n"); //img ref height
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int n0 = get_global_size(0);\n");
    source.append("     int n1 = get_global_size(1);\n");

    source.append("     " + input_type + " xt = xo[i + j*n0];\n");
    source.append("     " + input_type + " yt = yo[i + j*n0];\n");
    source.append("     int x = floor(xt);\n");
    source.append("     int y = floor(yt);\n");
    source.append("     if(x >= 0 && x < w && y >= 0 && y < h - 1)\n");
    source.append("     {\n");
    source.append("         " + input_type + "4 iv = {imgr[x+y*w], imgr[x+1+y*w], imgr[x+(y+1)*w], imgr[x+1+(y+1)*w]};\n");
    source.append("         " + input_type + " dx = xt - (" + input_type + ")x;\n");
    source.append("         " + input_type + " dy = yt - (" + input_type + ")y;\n");
    source.append("         " + input_type + " dxdy = dx*dy;\n");
    source.append("         " + input_type + " r = iv.s0*(1-dx-dy+dxdy) + iv.s1*(dx-dxdy) + iv.s2*(dy-dxdy) + iv.s3*dxdy;\n");
    source.append("         imgo[i + j*n0] = r;\n");
    source.append("     }\n");
    source.append("     else if(x >= 0 && x < w && y == h - 1)\n"); // border case
    source.append("     {\n");
    source.append("         " + input_type + "2 iv = {imgr[x+y*w], imgr[x+1+y*w]};\n");
    source.append("         " + input_type + " dx = xt - (" + input_type + ")x;\n");
    source.append("         " + input_type + " r = iv.s0*(1-dx) + iv.s1*(dx);\n");
    source.append("         imgo[i + j*n0] = r;\n");
    source.append("     };\n");

    // source.append("         " + input_type + "4 r = {imgr[x+y*w], imgr[x+1+y*w], imgr[x+(y+1)*w], imgr[x+1+(y+1)*w]};\n");
    // source.append("         imgo[i + j*n0] = bilinear(r,dx,dy);\n");
    source.append("};");
    return source;
};

std::string kernel_linear_interpolation_3d(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_linear_interpolation_3d(");
    source.append("     __global const " + input_type + " * xo,");
    source.append("     __global const " + input_type + " * yo,");
    source.append("     __global const " + input_type + " * zo,");
    source.append("     __global const " + input_type + " * imgr,");
    source.append("     __global " + input_type + " * imgo,");
    source.append("     int w,");   //img ref width
    source.append("     int h,");   //img ref height
    source.append("     int l)\n"); //img ref length
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int k = get_global_id(2);\n");
    source.append("     int n0 = get_global_size(0);\n");
    source.append("     int n1 = get_global_size(1);\n");
    source.append("     int n2 = get_global_size(2);\n");

    source.append("     " + input_type + " xt = xo[i + j*n0 + k*n0*n1];\n");
    source.append("     " + input_type + " yt = yo[i + j*n0 + k*n0*n1];\n");
    source.append("     " + input_type + " zt = zo[i + j*n0 + k*n0*n1];\n");
    source.append("     int x = floor(xt);\n");
    source.append("     int y = floor(yt);\n");
    source.append("     int z = floor(zt);\n");
    source.append("     if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l-1)\n");
    source.append("     {\n");
    source.append("         " + input_type + "4 iv = {imgr[x+y*w+z*w*h], imgr[x+1+y*w+z*w*h], imgr[x+(y+1)*w+z*w*h], imgr[x+1+(y+1)*w+z*w*h]};\n");
    source.append("         " + input_type + "4 iw = {imgr[x+y*w+(z+1)*w*h], imgr[x+1+y*w+(z+1)*w*h], imgr[x+(y+1)*w+(z+1)*w*h], imgr[x+1+(y+1)*w+(z+1)*w*h]};\n");
    source.append("         " + input_type + " dx = xt - (" + input_type + ")x;\n");
    source.append("         " + input_type + " dy = yt - (" + input_type + ")y;\n");
    source.append("         " + input_type + " dxdy = dx*dy;\n");
    source.append("         " + input_type + " rv = iv.s0*(1-dx-dy+dxdy) + iv.s1*(dx-dxdy) + iv.s2*(dy-dxdy) + iv.s3*dxdy;\n");
    source.append("         " + input_type + " rw = iw.s0*(1-dx-dy+dxdy) + iw.s1*(dx-dxdy) + iw.s2*(dy-dxdy) + iw.s3*dxdy;\n");
    source.append("         " + input_type + " dz = zt - (" + input_type + ")z;\n");
    source.append("         " + input_type + " r = rv*(1-dz) + rw*dz;\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = r;\n");
    source.append("     }\n");
    source.append("     else if(x >= 0 && x < w && y >= 0 && y < h && z == l-1)\n"); // border case
    source.append("     {\n");
    source.append("         " + input_type + "4 iv = {imgr[x+y*w+z*w*h], imgr[x+1+y*w+z*w*h], imgr[x+(y+1)*w+z*w*h], imgr[x+1+(y+1)*w+z*w*h]};\n");
    source.append("         " + input_type + " dx = xt - (" + input_type + ")x;\n");
    source.append("         " + input_type + " dy = yt - (" + input_type + ")y;\n");
    source.append("         " + input_type + " dxdy = dx*dy;\n");
    source.append("         " + input_type + " rv = iv.s0*(1-dx-dy+dxdy) + iv.s1*(dx-dxdy) + iv.s2*(dy-dxdy) + iv.s3*dxdy;\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = rv;\n");
    source.append("     };\n");
    source.append("};");
    return source;
};

std::string kernel_gradientx(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_gradientx(");
    source.append("     __global const " + input_type + " * imgr,");
    source.append("     __global " + input_type + " * imgo)\n");
    // source.append("     int w,");   //img ref width
    // source.append("     int h)\n"); //img ref height
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int k = get_global_id(2);\n");
    source.append("     int n0 = get_global_size(0);\n");
    source.append("     int n1 = get_global_size(1);\n");
    source.append("     int n2 = get_global_size(2);\n");

    source.append("     if(i == 0)\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = imgr[i+1 + j*n0 + k*n0*n1] - imgr[i + j*n0 + k*n0*n1];\n");
    source.append("     }\n");
    source.append("     else if(i == n0 - 1)\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = imgr[i + j*n0 + k*n0*n1] - imgr[i-1 + j*n0 + k*n0*n1];\n");
    source.append("     }\n");
    source.append("     else\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = 0.5*imgr[i+1 + j*n0 + k*n0*n1] - 0.5*imgr[i-1 + j*n0 + k*n0*n1];\n");
    source.append("     }\n");
    source.append("};");
    return source;
};

std::string kernel_gradienty(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_gradienty(");
    source.append("     __global const " + input_type + " * imgr,");
    source.append("     __global " + input_type + " * imgo)\n");
    // source.append("     int w,");   //img ref width
    // source.append("     int h)\n"); //img ref height
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int k = get_global_id(2);\n");
    source.append("     int n0 = get_global_size(0);\n");
    source.append("     int n1 = get_global_size(1);\n");
    source.append("     int n2 = get_global_size(2);\n");

    source.append("     if(j == 0)\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = imgr[i + (j+1)*n0 + k*n0*n1] - imgr[i + j*n0 + k*n0*n1];\n");
    source.append("     }\n");
    source.append("     else if(j == n1 - 1)\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = imgr[i + j*n0 + k*n0*n1] - imgr[i + (j-1)*n0 + k*n0*n1];\n");
    source.append("     }\n");
    source.append("     else\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = 0.5*imgr[i + (j+1)*n0 + k*n0*n1] - 0.5*imgr[i + (j-1)*n0 + k*n0*n1];\n");
    source.append("     }\n");
    source.append("};");
    return source;
};

std::string kernel_gradientz(std::string input_type)
{
    std::string source;
    source.append("__kernel void kernel_gradientz(");
    source.append("     __global const " + input_type + " * imgr,");
    source.append("     __global " + input_type + " * imgo)\n");
    // source.append("     int w,");   //img ref width
    // source.append("     int h)\n"); //img ref height
    source.append("{\n");
    source.append("     int i = get_global_id(0);\n");
    source.append("     int j = get_global_id(1);\n");
    source.append("     int k = get_global_id(2);\n");
    source.append("     int n0 = get_global_size(0);\n");
    source.append("     int n1 = get_global_size(1);\n");
    source.append("     int n2 = get_global_size(2);\n");

    source.append("     if(k == 0)\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = imgr[i + j*n0 + (k+1)*n0*n1] - imgr[i + j*n0 + k*n0*n1];\n");
    source.append("     }\n");
    source.append("     else if(k == n2 - 1)\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = imgr[i + j*n0 + k*n0*n1] - imgr[i + j*n0 + (k-1)*n0*n1];\n");
    source.append("     }\n");
    source.append("     else\n");
    source.append("     {\n");
    source.append("         imgo[i + j*n0 + k*n0*n1] = 0.5*imgr[i + j*n0 + (k+1)*n0*n1] - 0.5*imgr[i + j*n0 + (k-1)*n0*n1];\n");
    source.append("     }\n");
    source.append("};");
    return source;
};

}; //end namespace

#endif