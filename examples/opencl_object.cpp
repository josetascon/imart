/*
* @Author: Jose Tascon
* @Date:   2020-06-19 20:05:07
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-03-24 19:11:26
*/


// std libs
#include <iostream>

// local libs
#include "../src/opencl_object.h"
#include "../src/vector_opencl.h"
#include "../src/kernels.h"

using namespace imart;

int main()
{
    using type = int;
    // opencl_object cl_manager;

    // Manual set platform
    // cl_manager.set_platform_device(0,0);
    cl_manager.print_device_name();

    int N = 64;
    type scalar = 2;
    std::vector<type> vec_in(N);
    std::vector<type> vec_out(N);
    for(int i = 0; i < vec_in.size(); i++){ vec_in[i] = i; vec_out[i] = 1; };
    
    int err = 0;
    cl::Buffer in_buffer(cl_manager.get_context(), CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR, sizeof(type)*vec_in.size(), vec_in.data(), &err);
    std::cout << "[Status][in_buffer] Error: " << err << "\n";
    cl::Buffer out_buffer(cl_manager.get_context(), CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY, sizeof(type)*vec_out.size(), nullptr, &err);
    std::cout << "[Status][out_buffer] Error: " << err << "\n";

    std::string str_kernel = kernel_scalar( string_type<type>(), "+" );
    std::cout << str_kernel << std::endl;
    cl_manager.program(str_kernel, "kernel_scalar");
    cl_manager.arguments(in_buffer, out_buffer, scalar);
    cl_manager.execute(vec_in.size());
    
    cl::CommandQueue queue = cl_manager.get_queue();
    queue.enqueueReadBuffer(out_buffer, CL_TRUE, 0, sizeof(type)*vec_out.size(), vec_out.data());
    
    std::cout << "vec_in:\n[ ";
    for(int i = 0; i < vec_in.size(); i++){ std::cout << vec_in[i] << " "; };
    std::cout << "]" << std::endl;

    std::cout << "vec_out:\n[ ";
    for(int i = 0; i < vec_out.size(); i++){ std::cout << vec_out[i] << " "; };
    std::cout << "]" << std::endl;

    return 0;
}

    