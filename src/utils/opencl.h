/*
* @Author: jose
* @Date:   2020-01-27 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-01-27 00:00:00
*/

#ifndef __UTILS_OPENCL_H__
#define __UTILS_OPENCL_H__

// std libs
#include <iostream>     // std::cout
#include <string>       // std::string
#include <cassert>      // assert

// boost libs
#include <boost/stacktrace.hpp>

// opencl libs
#include <CL/cl.hpp>

namespace imart
{

// #ifndef NDEBUG
// #   define imart_assert_cl(Expr, Msg) \
//     imart_assert_opencl_error(#Expr, Expr, __FILE__, __LINE__, Msg)
// #else
// #   define imart_assert_cl(Expr, Msg) ;
// #endif

# define imart_assert_cl(status, msg) \
    imart_assert_opencl_error(status, __FILE__, __LINE__, msg)

std::string opencl_error(cl_int error);

void imart_assert_opencl_error(cl_int status, const char* file, int line, const char* msg)
{
    if (status != 0)
    {
        std::cerr << "\n******* OpenCL Error *******"
                  << "\n[Error] Information:\t" << msg
                  << "\n[Error] Error code:\t" << status
                  << "\n[Error] Description:\t" << opencl_error(status)
                  << "\n[Error] File:\t\t" << file
                  << "\n[Error] Line:\t\t" << line// << std::endl;
                  << "\n[Error] Backtrace:\n" << boost::stacktrace::stacktrace() << std::endl;
        assert(status == 0);
    };
};

// Print the error associciated with an error code
std::string opencl_error(cl_int error)
{
    // Print error message
    std::string error_code;
    switch(error)
    {
    case -1:
        error_code = "CL_DEVICE_NOT_FOUND";
        break;
    case -2:
        error_code = "CL_DEVICE_NOT_AVAILABLE";
        break;
    case -3:
        error_code = "CL_COMPILER_NOT_AVAILABLE";
        break;
    case -4:
        error_code = "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        break;
    case -5:
        error_code = "CL_OUT_OF_RESOURCES";
        break;
    case -6:
        error_code = "CL_OUT_OF_HOST_MEMORY";
        break;
    case -7:
        error_code = "CL_PROFILING_INFO_NOT_AVAILABLE";
        break;
    case -8:
        error_code = "CL_MEM_COPY_OVERLAP";
        break;
    case -9:
        error_code = "CL_IMAGE_FORMAT_MISMATCH";
        break;
    case -10:
        error_code = "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        break;
    case -11:
        error_code = "CL_BUILD_PROGRAM_FAILURE";
        break;
    case -12:
        error_code = "CL_MAP_FAILURE";
        break;
    case -13:
        error_code = "CL_MISALIGNED_SUB_BUFFER_OFFSET";
        break;
    case -14:
        error_code = "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        break;

    case -30:
        error_code = "CL_INVALID_VALUE";
        break;
    case -31:
        error_code = "CL_INVALID_DEVICE_TYPE";
        break;
    case -32:
        error_code = "CL_INVALID_PLATFORM";
        break;
    case -33:
        error_code = "CL_INVALID_DEVICE";
        break;
    case -34:
        error_code = "CL_INVALID_CONTEXT";
        break;
    case -35:
        error_code = "CL_INVALID_QUEUE_PROPERTIES";
        break;
    case -36:
        error_code = "CL_INVALID_COMMAND_QUEUE";
        break;
    case -37:
        error_code = "CL_INVALID_HOST_PTR";
        break;
    case -38:
        error_code = "CL_INVALID_MEM_OBJECT";
        break;
    case -39:
        error_code = "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        break;
    case -40:
        error_code = "CL_INVALID_IMAGE_SIZE";
        break;
    case -41:
        error_code = "CL_INVALID_SAMPLER";
        break;
    case -42:
        error_code = "CL_INVALID_BINARY";
        break;
    case -43:
        error_code = "CL_INVALID_BUILD_OPTIONS";
        break;
    case -44:
        error_code = "CL_INVALID_PROGRAM";
        break;
    case -45:
        error_code = "CL_INVALID_PROGRAM_EXECUTABLE";
        break;
    case -46:
        error_code = "CL_INVALID_KERNEL_NAME";
        break;
    case -47:
        error_code = "CL_INVALID_KERNEL_DEFINITION";
        break;
    case -48:
        error_code = "CL_INVALID_KERNEL";
        break;
    case -49:
        error_code = "CL_INVALID_ARG_INDEX";
        break;
    case -50:
        error_code = "CL_INVALID_ARG_VALUE";
        break;
    case -51:
        error_code = "CL_INVALID_ARG_SIZE";
        break;
    case -52:
        error_code = "CL_INVALID_KERNEL_ARGS";
        break;
    case -53:
        error_code = "CL_INVALID_WORK_DIMENSION";
        break;
    case -54:
        error_code = "CL_INVALID_WORK_GROUP_SIZE";
        break;
    case -55:
        error_code = "CL_INVALID_WORK_ITEM_SIZE";
        break;
    case -56:
        error_code = "CL_INVALID_GLOBAL_OFFSET";
        break;
    case -57:
        error_code = "CL_INVALID_EVENT_WAIT_LIST";
        break;
    case -58:
        error_code = "CL_INVALID_EVENT";
        break;
    case -59:
        error_code = "CL_INVALID_OPERATION";
        break;
    case -60:
        error_code = "CL_INVALID_GL_OBJECT";
        break;
    case -61:
        error_code = "CL_INVALID_BUFFER_SIZE";
        break;
    case -62:
        error_code = "CL_INVALID_MIP_LEVEL";
        break;
    case -63:
        error_code = "CL_INVALID_GLOBAL_WORK_SIZE";
        break;
    default:
        error_code = "UNRECOGNIZED ERROR CODE " + std::to_string(error);
    }
    return error_code;
}

}; //end namespace

#endif