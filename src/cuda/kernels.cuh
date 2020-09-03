/*
* @Author: jose
* @Date:   2020-08-24 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-08-24 00:00:00
*/

#ifndef __KERNELS_CUH__
#define __KERNELS_CUH__

// c lib
#include <stdio.h>

// std libs
#include <iostream>
#include <vector>
#include <cmath>

// local libs
// #include "interface.cuh"


// ===========================================
// Kernels
// ===========================================
template <typename type>
void cuda_kernel_assign( std::vector<int> & grid, std::vector<int> & block, 
                         type * vin, type value, int n );

template <typename type>
void cuda_kernel_copy( std::vector<int> & grid, std::vector<int> & block, 
                       const type * vin, type * vout, int n );

template <typename typein, typename typeout>
void cuda_kernel_cast( std::vector<int> & grid, std::vector<int> & block, 
                       const typein * vin, typeout * vout, int n );

// ===========================================
// Vector Kernels
// ===========================================
template <typename type>
void cuda_kernel_add_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n );

template <typename type>
void cuda_kernel_sub_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n );

template <typename type>
void cuda_kernel_sub_scalar_inv( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n );

template <typename type>
void cuda_kernel_mul_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n );

template <typename type>
void cuda_kernel_div_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n );

template <typename type>
void cuda_kernel_div_scalar_inv( std::vector<int> & grid, std::vector<int> & block, 
                                 const type * vin, type * vout, type scalar, int n );

template <typename type>
void cuda_kernel_pow_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n );

template <typename type>
void cuda_kernel_pow_scalar_inv( std::vector<int> & grid, std::vector<int> & block, 
                                 const type * vin, type * vout, type scalar, int n );

template <typename type>
void cuda_kernel_add( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n );

template <typename type>
void cuda_kernel_sub( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n );

template <typename type>
void cuda_kernel_mul( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n );

template <typename type>
void cuda_kernel_div( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n );

template <typename type>
void cuda_kernel_pow( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n );

// ===========================================
// Reduction Kernels
// ===========================================
template <typename type>
void cuda_kernel_sum( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin, type * vout, int n);

template <typename type>
void cuda_kernel_min( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin, type * vout, int n );

template <typename type>
void cuda_kernel_max( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin, type * vout, int n );

// ===========================================
// Image Kernels
// ===========================================
template <typename type>
void cuda_kernel_pad_2d( std::vector<int> & grid, std::vector<int> & block, 
                         const type * vin, type * vout, int start0, int start1,
                         int end0, int end1, int n0, int n1);

template <typename type>
void cuda_kernel_unpad_2d( std::vector<int> & grid, std::vector<int> & block, 
                            const type * vin, type * vout, int start0, int start1,
                            int end0, int end1, int n0, int n1);

template <typename type>
void cuda_kernel_pad_3d( std::vector<int> & grid, std::vector<int> & block, 
                         const type * vin, type * vout, int start0, int start1, int start2,
                         int end0, int end1, int end2, int n0, int n1, int n2 );

template <typename type>
void cuda_kernel_unpad_3d( std::vector<int> & grid, std::vector<int> & block, 
                            const type * vin, type * vout, int start0, int start1, int start2,
                            int end0, int end1, int end2, int n0, int n1, int n2 );

template <typename type>
void cuda_kernel_grid_2d( std::vector<int> & grid, std::vector<int> & block, 
                          type * x, type * y, double * sod, 
                          int n0, int n1);

template <typename type>
void cuda_kernel_grid_3d( std::vector<int> & grid, std::vector<int> & block,
                          type * x, type * y, type * z, double * sod, 
                          int n0, int n1, int n2 );

template <typename type>
void cuda_kernel_affine_2d( std::vector<int> & grid, std::vector<int> & block, 
                            const type * xin, const type * yin, 
                            type * xout, type * yout, 
                            const type * param, int n);

template <typename type>
void cuda_kernel_affine_3d( std::vector<int> & grid, std::vector<int> & block,
                            const type * xin, const type * yin, const type * zin,
                            type * xout, type * yout, type * zout,
                            const type * param, int n) ;

template <typename type>
void cuda_kernel_affine_sod_2d( std::vector<int> & grid, std::vector<int> & block,
                                const type * xin, const type * yin,
                                type * xout, type * yout,
                                const double * sod, int n);

template <typename type>
void cuda_kernel_affine_sod_3d( std::vector<int> & grid, std::vector<int> & block,
                                const type * xin, const type * yin, const type * zin,
                                type * xout, type * yout, type * zout,
                                const double * sod, int n);

template <typename type>
void cuda_kernel_dfield_2d( std::vector<int> & grid, std::vector<int> & block,
                            const type * xin, const type * yin,   // grid coordinates
                            const type * x, const type * y,       // vector field
                            type * xout, type * yout, int n);

template <typename type>
void cuda_kernel_dfield_3d( std::vector<int> & grid, std::vector<int> & block,
                            const type * xin, const type * yin, const type * zin, // grid coordinates
                            const type * x, const type * y, const type * z,       // vector field
                            type * xout, type * yout, type * zout,                // output coordinates
                            int n);

template <typename type>
void cuda_kernel_nearest_interpolation_2d( std::vector<int> & grid, std::vector<int> & block,
                                           const type * xo, const type * yo,
                                           const type * imgr, type * imgo,
                                           int w, int h,   //img ref width and height
                                           int n0, int n1); //img out dims

template <typename type>
void cuda_kernel_nearest_interpolation_3d( std::vector<int> & grid, std::vector<int> & block, 
                                           const type * xo, const type * yo, const type * zo, 
                                           const type * imgr, type * imgo,
                                           int w, int h, int l,    //img ref width, height and length
                                           int n0, int n1, int n2 );

template <typename type>
void cuda_kernel_linear_interpolation_2d( std::vector<int> & grid, std::vector<int> & block,
                                          const type * xo, const type * yo,
                                          const type * imgr, type * imgo,
                                          int w, int h,   //img ref width and height
                                          int n0, int n1); //img out dims

template <typename type>
void cuda_kernel_linear_interpolation_3d( std::vector<int> & grid, std::vector<int> & block,
                                          const type * xo, const type * yo, const type * zo,
                                          const type * imgr, type * imgo,
                                          int w, int h, int l, //img ref width, height and length
                                          int n0, int n1, int n2 );

template <typename type>
void cuda_kernel_gradientx( std::vector<int> & grid, std::vector<int> & block,
                            const type * imgr, type * imgo, 
                            int n0, int n1, int n2 );

template <typename type>
void cuda_kernel_gradienty( std::vector<int> & grid, std::vector<int> & block,
                            const type * imgr, type * imgo, 
                            int n0, int n1, int n2 );

template <typename type>
void cuda_kernel_gradientz( std::vector<int> & grid, std::vector<int> & block,
                            const type * imgr, type * imgo, 
                            int n0, int n1, int n2 );

template <typename type>
void cuda_kernel_convolution_2d( std::vector<int> & grid, std::vector<int> & block,
                                 const type * imgr, const type * kern, //kernel width
                                 type * imgo, int kwidth, int n0, int n1);

template <typename type>
void cuda_kernel_convolution_3d( std::vector<int> & grid, std::vector<int> & block,
                                 const type * imgr, const type * kern, //kernel width
                                 type * imgo, int kwidth, int n0, int n1, int n2 );

//NEW KERNELS
// template <typename type>
// void cuda_kernel_ ;

// template <typename type>
// void cuda_kernel_ ;

// template <typename type>
// void cuda_kernel_ ;

// template <typename type>
// void cuda_kernel_ ;

#endif