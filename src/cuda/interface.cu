/*
* @Author: jose
* @Date:   2020-08-24 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-08-24 00:00:00
*/

// local libs
#include "interface.cuh"

// ===========================================
// Check Errors
// ===========================================
#define imart_assert_cuda(status, msg) \
    imart_assert_cuda_error((status), __FILE__, __LINE__, msg);

void imart_assert_cuda_error(cudaError_t code, const char *file, int line, const char* msg, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"\n******* CUDA Error *******"\
                    "\n[Error] Information:\t%s"\
                    "\n[Error] Error code:\t%i"\
                    "\n[Error] Description:\t%s"\
                    "\n[Error] File:\t\t%s"\
                    "\n[Error] Line:\t\t%d\n",
                    msg, code, cudaGetErrorString(code), file, line);
      if (abort) exit(code); 
   };
};

// ===========================================
// Kernels
// ===========================================
// __global__ void kernel_print(const char * msg)
__global__ void kernel_print()
{
    printf("[GPU] Hello!\n");
    // printf("[GPU] %s\n", msg);
};

// ===========================================
// Functions
// ===========================================
void cuda_check_gpu()
{
    int devicesCount;
    cudaGetDeviceCount(&devicesCount);
    for(int deviceIndex = 0; deviceIndex < devicesCount; ++deviceIndex)
    {
        cudaDeviceProp deviceProperties;
        cudaGetDeviceProperties(&deviceProperties, deviceIndex);
        std::cout << "CUDA Device:\t" << deviceProperties.name << std::endl;
    }
};

void cuda_print()
{
    kernel_print<<<1, 1>>>();
    imart_assert_cuda( cudaPeekAtLastError(), "Fail to run kernel print" );
    imart_assert_cuda( cudaDeviceSynchronize(), "Fail to sync kernel print");
};

template <typename type>
void cuda_create_memory(type * & x, int size)
{
    imart_assert_cuda ( cudaMalloc(&x, size*sizeof(type)), "Memory allocation" ); 
};

template <typename type>
void cuda_clean_memory(type * & x)
{
    imart_assert_cuda( cudaFree(x), "Memory free" ); 
};

template <typename type>
void cuda_push_memory(type * x, type * data, int size, int offset)
{
    // printf("vector in:\n");
    // for(int i = 0; i < size; i++)
    //     printf("%f ",data[i]);
    imart_assert_cuda( cudaMemcpy(x, data, size*sizeof(type), cudaMemcpyHostToDevice), "Memory copy host to device" );
    // cudaMemcpy(x + offset, data, size*sizeof(type), cudaMemcpyHostToDevice);
};

template <typename type>
void cuda_push_memory(type * x, const type * data, int size, int offset)
{
    imart_assert_cuda( cudaMemcpy(x, data, size*sizeof(type), cudaMemcpyHostToDevice), "Memory copy host to device" );
};

template <typename type>
void cuda_pull_memory(type * x, type * data, int size, int offset)
{   
    // printf("pull\n");
    imart_assert_cuda( cudaMemcpy(data, x, size*sizeof(type), cudaMemcpyDeviceToHost), "Memory copy device to host" );
    // cudaMemcpy(data, x + offset, size*sizeof(type), cudaMemcpyDeviceToHost);
    // printf("vector out:\n");
    // for(int i = 0; i < size; i++)
    //     printf("%f ",data[i]);
};

// ===========================================
// Explicit instanciation
// ===========================================
template void cuda_create_memory<float>(float * & x, int size);
template void cuda_clean_memory<float>(float * & x);
template void cuda_push_memory<float>(float * x, float * data, int size, int offset);
template void cuda_push_memory<float>(float * x, const float * data, int size, int offset);
template void cuda_pull_memory<float>(float * x, float * data, int size, int offset);

template void cuda_create_memory<double>(double * & x, int size);
template void cuda_clean_memory<double>(double * & x);
template void cuda_push_memory<double>(double * x, double * data, int size, int offset);
template void cuda_push_memory<double>(double * x, const double * data, int size, int offset);
template void cuda_pull_memory<double>(double * x, double * data, int size, int offset);

template void cuda_create_memory<unsigned int>(unsigned int * & x, int size);
template void cuda_clean_memory<unsigned int>(unsigned int * & x);
template void cuda_push_memory<unsigned int>(unsigned int * x, unsigned int * data, int size, int offset);
template void cuda_push_memory<unsigned int>(unsigned int * x, const unsigned int * data, int size, int offset);
template void cuda_pull_memory<unsigned int>(unsigned int * x, unsigned int * data, int size, int offset);

template void cuda_create_memory<int>(int * & x, int size);
template void cuda_clean_memory<int>(int * & x);
template void cuda_push_memory<int>(int * x, int * data, int size, int offset);
template void cuda_push_memory<int>(int * x, const int * data, int size, int offset);
template void cuda_pull_memory<int>(int * x, int * data, int size, int offset);

template void cuda_create_memory<unsigned short>(unsigned short * & x, int size);
template void cuda_clean_memory<unsigned short>(unsigned short * & x);
template void cuda_push_memory<unsigned short>(unsigned short * x, unsigned short * data, int size, int offset);
template void cuda_push_memory<unsigned short>(unsigned short * x, const unsigned short * data, int size, int offset);
template void cuda_pull_memory<unsigned short>(unsigned short * x, unsigned short * data, int size, int offset);

template void cuda_create_memory<short>(short * & x, int size);
template void cuda_clean_memory<short>(short * & x);
template void cuda_push_memory<short>(short * x, short * data, int size, int offset);
template void cuda_push_memory<short>(short * x, const short * data, int size, int offset);
template void cuda_pull_memory<short>(short * x, short * data, int size, int offset);

template void cuda_create_memory<unsigned char>(unsigned char * & x, int size);
template void cuda_clean_memory<unsigned char>(unsigned char * & x);
template void cuda_push_memory<unsigned char>(unsigned char * x, unsigned char * data, int size, int offset);
template void cuda_push_memory<unsigned char>(unsigned char * x, const unsigned char * data, int size, int offset);
template void cuda_pull_memory<unsigned char>(unsigned char * x, unsigned char * data, int size, int offset);

template void cuda_create_memory<char>(char * & x, int size);
template void cuda_clean_memory<char>(char * & x);
template void cuda_push_memory<char>(char * x, char * data, int size, int offset);
template void cuda_push_memory<char>(char * x, const char * data, int size, int offset);
template void cuda_pull_memory<char>(char * x, char * data, int size, int offset);

