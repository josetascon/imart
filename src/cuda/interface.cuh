/*
* @Author: jose
* @Date:   2020-08-24 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-08-24 00:00:00
*/

#ifndef __INTERFACE_CUH__
#define __INTERFACE_CUH__

// c lib
#include <stdio.h>

// std libs
#include <iostream>
#include <vector>

// // ===========================================
// // Macro
// // ===========================================
// #define imart_assert_cuda(status, msg) \
//     imart_assert_opencl_error((status), __FILE__, __LINE__, msg);

// void imart_assert_opencl_error(cudaError_t code, const char *file, int line, const char* msg, bool abort=true);

// ===========================================
// Functions
// ===========================================
void cuda_check_gpu();

void cuda_print();

template <typename type>
void cuda_create_memory(type * & x, int size);

template <typename type>
void cuda_clean_memory(type * & x);

template <typename type>
void cuda_push_memory(type * x, type * ram, int size, int offset);

template <typename type>
void cuda_push_memory(type * x, const type * ram, int size, int offset);

template <typename type>
void cuda_pull_memory(type * x, type * ram, int size, int offset);

// template <typename type>
// void cuda_kernel_assign(std::vector<int> & grid, std::vector<int> & block, type * vin, type value, int n );

// template <typename type>
// void cuda_kernel_copy(std::vector<int> grid, std::vector<int> block, const type * vin, type * vout, int n );

#endif