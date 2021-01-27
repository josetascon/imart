/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __CONFIG_H__
#define __CONFIG_H__

// IMART precompiler definitions

#ifndef IMART_WITH_OPENMP
  #define IMART_WITH_OPENMP 1
#endif

#ifndef IMART_OPENMP_VECTOR_MIN_SIZE
  #define IMART_OPENMP_VECTOR_MIN_SIZE  5000
#endif

#ifndef IMART_WITH_OPENCL
  #define IMART_WITH_OPENCL 1
#endif

#ifndef IMART_WITH_CUDA
  #define IMART_WITH_CUDA 1
#endif

#ifndef IMART_WITH_FFTW
  #define IMART_WITH_FFTW 1
#endif

#ifndef IMART_WITH_CLFFT
  #define IMART_WITH_CLFFT 1
#endif

namespace imart
{
    
}; //end namespace

#endif