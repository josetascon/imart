/*
* @Author: jose
* @Date:   2020-08-24 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-08-24 00:00:00
*/

// local libs
#include "kernels.cuh"

// ===========================================
// Check Errors
// ===========================================
#define imart_assert_kernel(status, msg) \
    imart_assert_kernel_error((status), __FILE__, __LINE__, msg);

void imart_assert_kernel_error(cudaError_t code, const char *file, int line, const char* msg, bool abort=true)
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


// ===========================================
// Data Kernels
// ===========================================
template <typename type>
__global__ void kernel_assign(type * vin, type value, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vin[i] = value;
};

template <typename type>
__global__ void kernel_copy(const type * vin, type * vout, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = vin[i];
};

template <typename typein, typename typeout>
__global__ void kernel_cast(const typein * vin, typeout * vout, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = typeout(vin[i]);
};


// ===========================================
// Vector Kernels
// ===========================================
template <typename type>
__global__ void kernel_add_scalar(const type * vin, type * vout, type scalar, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = vin[i] + scalar;
};

template <typename type>
__global__ void kernel_sub_scalar(const type * vin, type * vout, type scalar, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = vin[i] - scalar;
};

template <typename type>
__global__ void kernel_sub_scalar_inv(const type * vin, type * vout, type scalar, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = scalar - vin[i];
};

template <typename type>
__global__ void kernel_mul_scalar(const type * vin, type * vout, type scalar, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = vin[i] * scalar;
};

template <typename type>
__global__ void kernel_div_scalar(const type * vin, type * vout, type scalar, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = vin[i] / scalar;
};

template <typename type>
__global__ void kernel_div_scalar_inv(const type * vin, type * vout, type scalar, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = scalar / vin[i];
};

template <typename type>
__global__ void kernel_pow_scalar(const type * vin, type * vout, type scalar, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = pow( vin[i], scalar );
};

template <typename type>
__global__ void kernel_pow_scalar_inv(const type * vin, type * vout, type scalar, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = pow( scalar, vin[i] );
};

template <typename type>
__global__ void kernel_add(const type * vin1, const type * vin2, type * vout, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = vin1[i] + vin2[i];
};

template <typename type>
__global__ void kernel_sub(const type * vin1, const type * vin2, type * vout, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = vin1[i] - vin2[i];
};

template <typename type>
__global__ void kernel_mul(const type * vin1, const type * vin2, type * vout, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = vin1[i] * vin2[i];
};

template <typename type>
__global__ void kernel_div(const type * vin1, const type * vin2, type * vout, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = vin1[i] / vin2[i];
};

template <typename type>
__global__ void kernel_pow(const type * vin1, const type * vin2, type * vout, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) vout[i] = pow( vin1[i], vin2[i] );
};

// ===========================================
// Reduction Kernels
// ===========================================
template <typename type>
__global__ void kernel_sum(const type *vin, type *vout, int n)
{
    __shared__ type sdata[256]; // Warning, threads should be 256
    unsigned int iii = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tid = threadIdx.x;
    type sum = 0;

    for (unsigned int i = iii; i < n; i += gridDim.x * blockDim.x)
    {
        sum += vin[i];
    };
    
    sdata[tid] = sum;
    __syncthreads();

    for (unsigned int s = blockDim.x >> 1; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            sdata[tid] += sdata[tid + s];
        };
        __syncthreads();
    };

    if (tid == 0) vout[blockIdx.x] = sdata[0];
};


template <typename type>
__global__ void kernel_min(const type *vin, type *vout, int n)
{
    __shared__ type sdata[256];
    unsigned int iii = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tid = threadIdx.x;
    type thread_result = vin[0];

    for (unsigned int i = iii; i < n; i += gridDim.x * blockDim.x)
    {
        type tmp = vin[i];
        thread_result = thread_result < tmp ? thread_result : tmp;
    };

    sdata[tid] = thread_result;
    __syncthreads();

    for (unsigned int s = blockDim.x >> 1; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            sdata[tid] =  sdata[tid] < sdata[tid + s]? sdata[tid] : sdata[tid + s];
        };
        __syncthreads();
    };

    if (tid == 0) vout[blockIdx.x] = sdata[0];
};

template <typename type>
__global__ void kernel_max(const type *vin, type *vout, int n)
{
    __shared__ type sdata[256];
    unsigned int iii = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tid = threadIdx.x;
    type thread_result = vin[0];

    for (unsigned int i = iii; i < n; i += gridDim.x * blockDim.x)
    {
        type tmp = vin[i];
        thread_result = thread_result > tmp ? thread_result : tmp;
    };

    sdata[tid] = thread_result;
    __syncthreads();

    for (unsigned int s = blockDim.x >> 1; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            sdata[tid] =  sdata[tid] > sdata[tid + s]? sdata[tid] : sdata[tid + s];
        };
        __syncthreads();
    };

    if (tid == 0) vout[blockIdx.x] = sdata[0];
};

// ===========================================
// Image Kernels
// ===========================================
template <typename type>
__global__ void kernel_pad_2d(const type * vin, type * vout, int start0, int start1,
                              int end0, int end1, int n0, int n1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;

    int wo = n0+start0+end0;

    if (i < n0 && j < n1) // width = n0, heigth = n1
    {
        vout[start0+i + (start1+j)*wo] = vin[i + j*n0];
    };
};

template <typename type>
__global__ void kernel_unpad_2d(const type * vin, type * vout, int start0, int start1,
                              int end0, int end1, int n0, int n1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;

    int wo = n0+start0+end0;

    if (i < n0 && j < n1) // width = n0, heigth = n1
    {
        vout[i + j*n0] = vin[start0+i + (start1+j)*wo];
    };
};

template <typename type>
__global__ void kernel_pad_3d(const type * vin, type * vout, int start0, int start1, int start2,
                              int end0, int end1, int end2, int n0, int n1, int n2)
{
    // int blockIdx_z = __float2int_rd(blockIdx.y * invBlocksInY);
    // int blockIdx_y = blockIdx.y - (blockIdx_z * blocksInY);
    // int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    // int j = (blockIdx_y * blockDim.y) + threadIdx.y;
    // int k = (blockIdx_z * blockDim.z) + threadIdx.z;
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y) + threadIdx.y;
    int k = (blockIdx.z * blockDim.z) + threadIdx.z;

    int wo = n0+start0+end0; //vout size
    int ho = n1+start1+end1; //vout size

    if (i < n0 && j < n1 && k < n2) // width = n0, height = n1, depth = n2
    {
        vout[start0+i + (start1+j)*wo + (start2+k)*wo*ho] = vin[i + j*n0 + k*n0*n1];
    };
};

template <typename type>
__global__ void kernel_unpad_3d(const type * vin, type * vout, int start0, int start1, int start2,
                              int end0, int end1, int end2, int n0, int n1, int n2)
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y) + threadIdx.y;
    int k = (blockIdx.z * blockDim.z) + threadIdx.z;

    int wo = n0+start0+end0; //vout size
    int ho = n1+start1+end1; //vout size

    if (i < n0 && j < n1 && k < n2) // width = n0, height = n1, depth = n2
    {
        vout[i + j*n0 + k*n0*n1] = vin[start0+i + (start1+j)*wo + (start2+k)*wo*ho];
    };
};

template <typename type>
__global__ void kernel_grid_2d( type * x, type * y, double * sod, 
                                int n0, int n1)
{
    // consider sod conversion to float to support all gpu
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;

    double c0 = sod[0]; double c1 = sod[1];
    double o0 = sod[2]; double o1 = sod[3];
    double d0 = sod[4]; double d1 = sod[5];
    double d2 = sod[6]; double d3 = sod[7];

    if (i < n0 && j < n1) // width = n0, heigth = n1
    {
        x[i+j*n0] = (type)(d0*c0*i + d1*c1*j + o0);
        y[i+j*n0] = (type)(d2*c0*i + d3*c1*j + o1);
    };
};

template <typename type>
__global__ void kernel_grid_3d( type * x, type * y, type * z, double * sod, 
                                int n0, int n1, int n2)
{
    // consider sod conversion to float to support all gpu
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y) + threadIdx.y;
    int k = (blockIdx.z * blockDim.z) + threadIdx.z;

    double c0 = sod[0]; double c1 = sod[1]; double c2 = sod[2];
    double o0 = sod[3]; double o1 = sod[4]; double o2 = sod[5];
    double d0 = sod[6]; double d1 = sod[7]; double d2 = sod[8];
    double d3 = sod[9]; double d4 = sod[10]; double d5 = sod[11];
    double d6 = sod[12]; double d7 = sod[13]; double d8 = sod[14];

    if (i < n0 && j < n1 && k < n2) // width = n0, height = n1, depth = n2
    {
        x[i + j*n0 + k*n0*n1] = (type)(d0*c0*i + d1*c1*j + d2*c2*k + o0);
        y[i + j*n0 + k*n0*n1] = (type)(d3*c0*i + d4*c1*j + d5*c2*k + o1);
        z[i + j*n0 + k*n0*n1] = (type)(d6*c0*i + d7*c1*j + d8*c2*k + o2);
    };
};

template <typename type>
__global__ void kernel_affine_2d( const type * xin, const type * yin, 
                                  type * xout, type * yout, 
                                  const type * param, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x; // one dimension, buffer in and out xy equal size
    type a0 = param[0]; type a1 = param[1];
    type a2 = param[2]; type a3 = param[3];
    type t0 = param[4]; type t1 = param[5];
    if (i < n)
    {
        xout[i] = (type)(a0*xin[i] + a1*yin[i] + t0);
        yout[i] = (type)(a2*xin[i] + a3*yin[i] + t1);
    };
};

template <typename type>
__global__ void kernel_affine_3d( const type * xin, const type * yin, const type * zin,
                                  type * xout, type * yout, type * zout,
                                  const type * param, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x; // one dimension, buffer in and out xyz equal size
    type a0 = param[0]; type a1 = param[1]; type a2 = param[2];
    type a3 = param[3]; type a4 = param[4]; type a5 = param[5];
    type a6 = param[6]; type a7 = param[7]; type a8 = param[8];
    type t0 = param[9]; type t1 = param[10]; type t2 = param[11];
    if (i < n)
    {
        xout[i] = (type)(a0*xin[i] + a1*yin[i] + a2*zin[i] + t0);
        yout[i] = (type)(a3*xin[i] + a4*yin[i] + a5*zin[i] + t1);
        zout[i] = (type)(a6*xin[i] + a7*yin[i] + a8*zin[i] + t2);
    };
};

template <typename type>
__global__ void kernel_affine_sod_2d( const type * xin, const type * yin,
                                      type * xout, type * yout,
                                      const double * sod, int n)
{
    // consider sod conversion to float to support all gpu
    int i = blockDim.x * blockIdx.x + threadIdx.x; // one dimension, buffer in and out xy equal size
    double c0 = sod[0]; double c1 = sod[1];
    double o0 = sod[2]; double o1 = sod[3];
    double d0 = sod[4]; double d1 = sod[5];
    double d2 = sod[6]; double d3 = sod[7];
    if (i < n)
    {
        xout[i] = (type)(d0*c0*xin[i] + d1*c1*yin[i] + o0);
        yout[i] = (type)(d2*c0*xin[i] + d3*c1*yin[i] + o1);
    }
};

template <typename type>
__global__ void kernel_affine_sod_3d( const type * xin, const type * yin, const type * zin,
                                      type * xout, type * yout, type * zout,
                                      const double * sod, int n)
{
    // consider sod conversion to float to support all gpu
    int i = blockDim.x * blockIdx.x + threadIdx.x; // one dimension, buffer in and out xyz equal size
    double c0 = sod[0]; double c1 = sod[1]; double c2 = sod[2];
    double o0 = sod[3]; double o1 = sod[4]; double o2 = sod[5];
    double d0 = sod[6]; double d1 = sod[7]; double d2 = sod[8];
    double d3 = sod[9]; double d4 = sod[10]; double d5 = sod[11];
    double d6 = sod[12]; double d7 = sod[13]; double d8 = sod[14];
    if (i < n)
    {
        xout[i] = (type)(d0*c0*xin[i] + d1*c1*yin[i] + d2*c2*zin[i] + o0);
        yout[i] = (type)(d3*c0*xin[i] + d4*c1*yin[i] + d5*c2*zin[i] + o1);
        zout[i] = (type)(d6*c0*xin[i] + d7*c1*yin[i] + d8*c2*zin[i] + o2);
    };
};

template <typename type>
__global__ void kernel_dfield_2d( const type * xin, const type * yin,   // grid coordinates
                                  const type * x, const type * y,       // vector field
                                  type * xout, type * yout, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x; // one dimension, buffer in and out xy equal size
    if (i < n)
    {
        xout[i] = xin[i] + x[i];
        yout[i] = yin[i] + y[i];
    };
};

template <typename type>
__global__ void kernel_dfield_3d( const type * xin, const type * yin, const type * zin, // grid coordinates
                                  const type * x, const type * y, const type * z,       // vector field
                                  type * xout, type * yout, type * zout,                // output coordinates
                                  int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x; // one dimension, buffer in and out xy equal size
    if (i < n)
    {
        xout[i] = xin[i] + x[i];
        yout[i] = yin[i] + y[i];
        zout[i] = zin[i] + z[i];
    };
};

template <typename type>
__global__ void kernel_nearest_interpolation_2d( const type * xo, const type * yo,
                                                 const type * imgr, type * imgo,
                                                 int w, int h,   //img ref width and height
                                                 int n0, int n1) //img out dims
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;

    if (i < n0 && j < n1)
    {
        int x = round(xo[i + j*n0]);
        int y = round(yo[i + j*n0]);
        if(x >= 0 && x < w && y >= 0 && y < h)
        {
            imgo[i + j*n0] = imgr[x + y*w];
        };
    };
};

template <typename type>
__global__ void kernel_nearest_interpolation_3d( const type * xo, const type * yo, const type * zo, 
                                                 const type * imgr, type * imgo,
                                                int w, int h, int l,    //img ref width, height and length
                                                int n0, int n1, int n2)
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y) + threadIdx.y;
    int k = (blockIdx.z * blockDim.z) + threadIdx.z;

    if (i < n0 && j < n1 && k < n2)
    {
        int x = round(xo[i + j*n0 + k*n0*n1]);
        int y = round(yo[i + j*n0 + k*n0*n1]);
        int z = round(yo[i + j*n0 + k*n0*n1]);
        if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l)
        {
            imgo[i + j*n0 + k*n0*n1] = imgr[x + y*w + z*w*h];
        };
    };
};

template <typename type>
__global__ void kernel_linear_interpolation_2d( const type * xo, const type * yo,
                                                const type * imgr, type * imgo,
                                                int w, int h,   //img ref width and height
                                                int n0, int n1) //img out dims
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;

    if (i < n0 && j < n1)
    {
        type xt = xo[i + j*n0];
        type yt = yo[i + j*n0];
        int x = floor(xt);
        int y = floor(yt);
        if(x >= 0 && x < w && y >= 0 && y < h - 1)
        {
            // __shared__ iv[4];
            type iv[4] = {imgr[x+y*w], imgr[x+1+y*w], imgr[x+(y+1)*w], imgr[x+1+(y+1)*w]};
            type dx = xt - (type)x;
            type dy = yt - (type)y;
            type dxdy = dx*dy;
            type r = iv[0]*(1-dx-dy+dxdy) + iv[1]*(dx-dxdy) + iv[2]*(dy-dxdy) + iv[3]*dxdy;
            imgo[i + j*n0] = r;
        }
        else if(x >= 0 && x < w && y == h - 1) // border case
        {
            type iv[2] = {imgr[x+y*w], imgr[x+1+y*w]};
            type dx = xt - (type)x;
            type r = iv[0]*(1-dx) + iv[1]*(dx);
            imgo[i + j*n0] = r;
        };
    };
};

template <typename type>
__global__ void kernel_linear_interpolation_3d( const type * xo, const type * yo, const type * zo,
                                                const type * imgr, type * imgo,
                                                int w, int h, int l, //img ref width, height and length
                                                int n0, int n1, int n2)
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y) + threadIdx.y;
    int k = (blockIdx.z * blockDim.z) + threadIdx.z;

    if (i < n0 && j < n1 && k < n2)
    {
        type xt = xo[i + j*n0 + k*n0*n1];
        type yt = yo[i + j*n0 + k*n0*n1];
        type zt = zo[i + j*n0 + k*n0*n1];
        int x = floor(xt);
        int y = floor(yt);
        int z = floor(zt);
        if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < l-1)
        {
            type iv[4] = {imgr[x+y*w+z*w*h], imgr[x+1+y*w+z*w*h], imgr[x+(y+1)*w+z*w*h], imgr[x+1+(y+1)*w+z*w*h]};
            type iw[4] = {imgr[x+y*w+(z+1)*w*h], imgr[x+1+y*w+(z+1)*w*h], imgr[x+(y+1)*w+(z+1)*w*h], imgr[x+1+(y+1)*w+(z+1)*w*h]};
            type dx = xt - (type)x;
            type dy = yt - (type)y;
            type dxdy = dx*dy;
            type rv = iv[0]*(1-dx-dy+dxdy) + iv[1]*(dx-dxdy) + iv[2]*(dy-dxdy) + iv[3]*dxdy;
            type rw = iw[0]*(1-dx-dy+dxdy) + iw[1]*(dx-dxdy) + iw[2]*(dy-dxdy) + iw[3]*dxdy;
            type dz = zt - (type)z;
            type r = rv*(1-dz) + rw*dz;
            imgo[i + j*n0 + k*n0*n1] = r;
        }
        else if(x >= 0 && x < w && y >= 0 && y < h && z == l-1) // border case
        {
            type iv[4] = {imgr[x+y*w+z*w*h], imgr[x+1+y*w+z*w*h], imgr[x+(y+1)*w+z*w*h], imgr[x+1+(y+1)*w+z*w*h]};
            type dx = xt - (type)x;
            type dy = yt - (type)y;
            type dxdy = dx*dy;
            type rv = iv[0]*(1-dx-dy+dxdy) + iv[1]*(dx-dxdy) + iv[2]*(dy-dxdy) + iv[3]*dxdy;
            imgo[i + j*n0 + k*n0*n1] = rv;
        };
    };
};

template <typename type>
__global__ void kernel_gradientx( const type * imgr, type * imgo, 
                                  int n0, int n1, int n2)
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y) + threadIdx.y;
    int k = (blockIdx.z * blockDim.z) + threadIdx.z;

    if (i < n0 && j < n1 && (k == 0 || k < n2))
    {
        if(i == 0)
        {
            imgo[i + j*n0 + k*n0*n1] = imgr[i+1 + j*n0 + k*n0*n1] - imgr[i + j*n0 + k*n0*n1];
        }
        else if(i == n0 - 1)
        {
            imgo[i + j*n0 + k*n0*n1] = imgr[i + j*n0 + k*n0*n1] - imgr[i-1 + j*n0 + k*n0*n1];
        }
        else
        {
            imgo[i + j*n0 + k*n0*n1] = 0.5*imgr[i+1 + j*n0 + k*n0*n1] - 0.5*imgr[i-1 + j*n0 + k*n0*n1];
        };
    };
};

template <typename type>
__global__ void kernel_gradienty( const type * imgr, type * imgo, 
                                  int n0, int n1, int n2)
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y) + threadIdx.y;
    int k = (blockIdx.z * blockDim.z) + threadIdx.z;

    if (i < n0 && j < n1 && (k == 0 || k < n2))
    {
        if(j == 0)
        {
            imgo[i + j*n0 + k*n0*n1] = imgr[i + (j+1)*n0 + k*n0*n1] - imgr[i + j*n0 + k*n0*n1];
        }
        else if(j == n1 - 1)
        {
            imgo[i + j*n0 + k*n0*n1] = imgr[i + j*n0 + k*n0*n1] - imgr[i + (j-1)*n0 + k*n0*n1];
        }
        else
        {
            imgo[i + j*n0 + k*n0*n1] = 0.5*imgr[i + (j+1)*n0 + k*n0*n1] - 0.5*imgr[i + (j-1)*n0 + k*n0*n1];
        };
    };
};

template <typename type>
__global__ void kernel_gradientz( const type * imgr, type * imgo, 
                                  int n0, int n1, int n2)
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y) + threadIdx.y;
    int k = (blockIdx.z * blockDim.z) + threadIdx.z;

    if (i < n0 && j < n1 && k < n2)
    {
        if(k == 0)
        {
            imgo[i + j*n0 + k*n0*n1] = imgr[i + j*n0 + (k+1)*n0*n1] - imgr[i + j*n0 + k*n0*n1];
        }
        else if(k == n2 - 1)
        {
            imgo[i + j*n0 + k*n0*n1] = imgr[i + j*n0 + k*n0*n1] - imgr[i + j*n0 + (k-1)*n0*n1];
        }
        else
        {
            imgo[i + j*n0 + k*n0*n1] = 0.5*imgr[i + j*n0 + (k+1)*n0*n1] - 0.5*imgr[i + j*n0 + (k-1)*n0*n1];
        };
    };
};

template <typename type>
__global__ void kernel_convolution_2d( const type * imgr, const type * kern, //kernel width
                                       type * imgo, int kwidth, int n0, int n1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;
    int off = floor(kwidth/2.0);

    if(i >= off && i < n0 - off && j >= off && j < n1 - off)
    {
        type sum = 0;
        for (int p = 0; p < kwidth; p++)
        {
            for (int q = 0; q < kwidth; q++)
            {
                sum += imgr[i+p-off + (j+q-off)*n0] * kern[p*kwidth + q];
            };
        };
        imgo[i + j*n0] = sum;
    };
};

template <typename type>
__global__ void kernel_convolution_3d( const type * imgr, const type * kern, //kernel width
                                       type * imgo, int kwidth, int n0, int n1, int n2)
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int j = (blockIdx.y * blockDim.y) + threadIdx.y;
    int k = (blockIdx.z * blockDim.z) + threadIdx.z;
    
    int off = floor(kwidth/2.0);

    if(i >= off && i < n0 - off && j >= off && j < n1 - off && k >= off && k < n2 - off)
    {
        type sum = 0;
        for (int r = 0; r < kwidth; r++)
        {
            for (int p = 0; p < kwidth; p++)
            {
                for (int q = 0; q < kwidth; q++)
                {
                    sum += imgr[i+p-off + (j+q-off)*n0 + (k-off)*n0*n1] * kern[r*kwidth*kwidth + p*kwidth + q];
                };
            };
        };
        imgo[i + j*n0 + k*n0*n1] = sum;
    };
};


















// ===========================================
// Kernels Calls
// ===========================================


// ===========================================
// Data Kernels
// ===========================================
template <typename type>
void cuda_kernel_assign( std::vector<int> & grid, std::vector<int> & block, 
                         type * vin, type value, int n )
{
    // printf("kernel assign init\n");
    // printf("block: [%i, %i, %i]\n", block[0], block[1] , block[2]);
    // printf("grid: [%i, %i, %i]\n", grid[0], grid[1] , grid[2]);
    // printf("address: %x\n", vin);
    // printf("value: %f\n", value);
    // printf("size: %i\n", n);

    dim3 grd(grid[0]);
    dim3 blk(block[0]);

    kernel_assign<<<grd,blk>>>(vin, value, n);
    // kernel_assign<type><<<grd,blk>>>(vin, value, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel assign" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel assign" );
    // printf("kernel assign finish\n");
};

template <typename type>
void cuda_kernel_copy( std::vector<int> & grid, std::vector<int> & block,
                       const type * vin, type * vout, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_copy<<<grd,blk>>>(vin, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel copy" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel copy" );
};

template <typename typein, typename typeout>
void cuda_kernel_cast( std::vector<int> & grid, std::vector<int> & block, 
                       const typein * vin, typeout * vout, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_cast<typein,typeout><<<grd,blk>>>(vin, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel cast" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel cast" );
};

// ===========================================
// Vector Kernels
// ===========================================
template <typename type>
void cuda_kernel_add_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_add_scalar<<<grd,blk>>>(vin, vout, scalar, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel add scalar" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel add scalar" );
};

template <typename type>
void cuda_kernel_sub_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_sub_scalar<<<grd,blk>>>(vin, vout, scalar, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel sub scalar" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel sub scalar" );
};

template <typename type>
void cuda_kernel_sub_scalar_inv( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_sub_scalar_inv<<<grd,blk>>>(vin, vout, scalar, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel sub scalar inv" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel sub scalar inv" );
};

template <typename type>
void cuda_kernel_mul_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_mul_scalar<<<grd,blk>>>(vin, vout, scalar, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel mul scalar" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel mul scalar" );
};

template <typename type>
void cuda_kernel_div_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_div_scalar<<<grd,blk>>>(vin, vout, scalar, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel div scalar" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel div scalar" );
};

template <typename type>
void cuda_kernel_div_scalar_inv( std::vector<int> & grid, std::vector<int> & block, 
                                 const type * vin, type * vout, type scalar, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_div_scalar_inv<<<grd,blk>>>(vin, vout, scalar, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel div scalar inv" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel div scalar inv" );
};

template <typename type>
void cuda_kernel_pow_scalar( std::vector<int> & grid, std::vector<int> & block, 
                             const type * vin, type * vout, type scalar, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_pow_scalar<<<grd,blk>>>(vin, vout, scalar, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel pow scalar" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel pow scalar" );
};

template <typename type>
void cuda_kernel_pow_scalar_inv( std::vector<int> & grid, std::vector<int> & block, 
                                 const type * vin, type * vout, type scalar, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_pow_scalar_inv<<<grd,blk>>>(vin, vout, scalar, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel pow scalar inv" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel pow scalar inv" );
};

template <typename type>
void cuda_kernel_add( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_add<<<grd,blk>>>(vin1, vin2, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel add" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel add" );
};

template <typename type>
void cuda_kernel_sub( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_sub<<<grd,blk>>>(vin1, vin2, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel sub" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel sub" );
};

template <typename type>
void cuda_kernel_mul( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_mul<<<grd,blk>>>(vin1, vin2, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel mul" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel mul" );
};

template <typename type>
void cuda_kernel_div( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_div<<<grd,blk>>>(vin1, vin2, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel div" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel div" );
};

template <typename type>
void cuda_kernel_pow( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin1, const type * vin2, type * vout, int n )
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_pow<<<grd,blk>>>(vin1, vin2, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel" );
};


// ===========================================
// Reduction Kernels
// ===========================================
template <typename type>
void cuda_kernel_sum( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin, type * vout, int n)
{
    // printf("kernel sum init\n");
    // printf("block: [%i, %i, %i]\n", block[0], block[1] , block[2]);
    // printf("grid: [%i, %i, %i]\n", grid[0], grid[1] , grid[2]);
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_sum<<<grd,blk>>>(vin, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel sum" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel sum" );
};

template <typename type>
void cuda_kernel_min( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin, type * vout, int n)
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_min<<<grd,blk>>>(vin, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel min" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel min" );
};

template <typename type>
void cuda_kernel_max( std::vector<int> & grid, std::vector<int> & block, 
                      const type * vin, type * vout, int n)
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_max<<<grd,blk>>>(vin, vout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel max" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel max" );
};

// ===========================================
// Image Kernels
// ===========================================
template <typename type>
void cuda_kernel_pad_2d( std::vector<int> & grid, std::vector<int> & block, 
                         const type * vin, type * vout, int start0, int start1,
                         int end0, int end1, int n0, int n1 )
{
    // printf("kernel pad init\n");
    // printf("block: [%i, %i, %i]\n", block[0], block[1] , block[2]);
    // printf("grid: [%i, %i, %i]\n", grid[0], grid[1] , grid[2]);
    dim3 grd(grid[0],grid[1]);
    dim3 blk(block[0],block[1]);
    kernel_pad_2d<<<grd,blk>>>(vin, vout, start0, start1, end0, end1, n0, n1);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel pad 2d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel pad 2d" );
};

template <typename type>
void cuda_kernel_unpad_2d( std::vector<int> & grid, std::vector<int> & block, 
                           const type * vin, type * vout, int start0, int start1,
                           int end0, int end1, int n0, int n1 )
{
    dim3 grd(grid[0],grid[1]);
    dim3 blk(block[0],block[1]);
    kernel_unpad_2d<<<grd,blk>>>(vin, vout, start0, start1, end0, end1, n0, n1);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel unpad 2d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel unpad 2d" );
};

template <typename type>
void cuda_kernel_pad_3d( std::vector<int> & grid, std::vector<int> & block, 
                         const type * vin, type * vout, int start0, int start1, int start2,
                         int end0, int end1, int end2, int n0, int n1, int n2 )
{
    // printf("kernel pad 3d init\n");
    // printf("block: [%i, %i, %i]\n", block[0], block[1] , block[2]);
    // printf("grid: [%i, %i, %i]\n", grid[0], grid[1] , grid[2]);
    dim3 grd(grid[0],grid[1],grid[2]);
    dim3 blk(block[0],block[1],block[2]);
    kernel_pad_3d<<<grd,blk>>>(vin, vout, start0, start1, start2, end0, end1, end2, n0, n1, n2);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel pad 3d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel pad 3d" );
};

template <typename type>
void cuda_kernel_unpad_3d( std::vector<int> & grid, std::vector<int> & block, 
                           const type * vin, type * vout, int start0, int start1, int start2,
                         int end0, int end1, int end2, int n0, int n1, int n2 )
{
    dim3 grd(grid[0],grid[1],grid[2]);
    dim3 blk(block[0],block[1],block[2]);
    kernel_unpad_3d<<<grd,blk>>>(vin, vout, start0, start1, start2, end0, end1, end2, n0, n1, n2);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel unpad 3d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel unpad 3d" );
};

template <typename type>
void cuda_kernel_grid_2d( std::vector<int> & grid, std::vector<int> & block, 
                          type * x, type * y, double * sod, 
                          int n0, int n1)
{
    dim3 grd(grid[0],grid[1]);
    dim3 blk(block[0],block[1]);
    kernel_grid_2d<<<grd,blk>>>(x, y, sod, n0, n1);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel grid 2d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel grid 2d" );
};

template <typename type>
void cuda_kernel_grid_3d( std::vector<int> & grid, std::vector<int> & block,
                          type * x, type * y, type * z, double * sod, 
                          int n0, int n1, int n2)
{
    dim3 grd(grid[0],grid[1],grid[2]);
    dim3 blk(block[0],block[1],block[2]);
    kernel_grid_3d<<<grd,blk>>>(x, y, z, sod, n0, n1, n2);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel grid 3d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel grid 3d" );
};

template <typename type>
void cuda_kernel_affine_2d( std::vector<int> & grid, std::vector<int> & block, 
                            const type * xin, const type * yin, 
                            type * xout, type * yout, 
                            const type * param, int n)
{
    // printf("kernel affine 2d init\n");
    // printf("block: [%i, %i, %i]\n", block[0], block[1] , block[2]);
    // printf("grid: [%i, %i, %i]\n", grid[0], grid[1] , grid[2]);
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_affine_2d<<<grd,blk>>>(xin, yin, xout, yout, param, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel affine 2d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel affine 2d" );
};

template <typename type>
void cuda_kernel_affine_3d( std::vector<int> & grid, std::vector<int> & block,
                            const type * xin, const type * yin, const type * zin,
                            type * xout, type * yout, type * zout,
                            const type * param, int n)
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_affine_3d<<<grd,blk>>>(xin, yin, zin, xout, yout, zout, param, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel affine 3d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel affine 3d" );
};

template <typename type>
void cuda_kernel_affine_sod_2d( std::vector<int> & grid, std::vector<int> & block,
                                const type * xin, const type * yin,
                                type * xout, type * yout,
                                const double * sod, int n)
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_affine_sod_2d<<<grd,blk>>>(xin, yin, xout, yout, sod, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel affine sod 2d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel affine sod 2d" );
};

template <typename type>
void cuda_kernel_affine_sod_3d( std::vector<int> & grid, std::vector<int> & block,
                                const type * xin, const type * yin, const type * zin,
                                type * xout, type * yout, type * zout,
                                const double * sod, int n)
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_affine_sod_3d<<<grd,blk>>>(xin, yin, zin, xout, yout, zout, sod, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel affine sod 3d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel affine sod 3d" );
};

template <typename type>
void cuda_kernel_dfield_2d( std::vector<int> & grid, std::vector<int> & block,
                            const type * xin, const type * yin,   // grid coordinates
                            const type * x, const type * y,       // vector field
                            type * xout, type * yout, int n)
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_dfield_2d<<<grd,blk>>>(xin, yin, x, y, xout, yout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel dfield 2d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel dfield 2d" );
};

template <typename type>
void cuda_kernel_dfield_3d( std::vector<int> & grid, std::vector<int> & block,
                            const type * xin, const type * yin, const type * zin, // grid coordinates
                            const type * x, const type * y, const type * z,       // vector field
                            type * xout, type * yout, type * zout,                // output coordinates
                            int n)
{
    dim3 grd(grid[0]);
    dim3 blk(block[0]);
    kernel_dfield_3d<<<grd,blk>>>(xin, yin, zin, x, y, z, xout, yout, zout, n);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel dfield 3d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel dfield 3d" );
};

template <typename type>
void cuda_kernel_nearest_interpolation_2d( std::vector<int> & grid, std::vector<int> & block,
                                           const type * xo, const type * yo,
                                           const type * imgr, type * imgo,
                                           int w, int h,   //img ref width and height
                                           int n0, int n1) //img out dims
{
    dim3 grd(grid[0],grid[1]);
    dim3 blk(block[0],block[1]);
    kernel_nearest_interpolation_2d<<<grd,blk>>>(xo, yo, imgr, imgo, w, h, n0, n1);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel nearest interpolation 2d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel nearest interpolation 2d" );
};

template <typename type>
void cuda_kernel_nearest_interpolation_3d( std::vector<int> & grid, std::vector<int> & block, 
                                           const type * xo, const type * yo, const type * zo, 
                                           const type * imgr, type * imgo,
                                           int w, int h, int l,    //img ref width, height and length
                                           int n0, int n1, int n2)
{
    dim3 grd(grid[0],grid[1],grid[2]);
    dim3 blk(block[0],block[1],block[2]);
    kernel_nearest_interpolation_3d<<<grd,blk>>>(xo, yo, zo, imgr, imgo, w, h, l, n0, n1, n2);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel nearest interpolation 3d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel nearest interpolation 3d" );
};

template <typename type>
void cuda_kernel_linear_interpolation_2d( std::vector<int> & grid, std::vector<int> & block,
                                          const type * xo, const type * yo,
                                          const type * imgr, type * imgo,
                                          int w, int h,   //img ref width and height
                                          int n0, int n1) //img out dims
{
    dim3 grd(grid[0],grid[1]);
    dim3 blk(block[0],block[1]);
    kernel_linear_interpolation_2d<<<grd,blk>>>(xo, yo, imgr, imgo, w, h, n0, n1);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel linear interpolation 2d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel linear interpolation 2d" );
};

template <typename type>
void cuda_kernel_linear_interpolation_3d( std::vector<int> & grid, std::vector<int> & block,
                                          const type * xo, const type * yo, const type * zo,
                                          const type * imgr, type * imgo,
                                          int w, int h, int l, //img ref width, height and length
                                          int n0, int n1, int n2)
{
    dim3 grd(grid[0],grid[1],grid[2]);
    dim3 blk(block[0],block[1],block[2]);
    kernel_linear_interpolation_3d<<<grd,blk>>>(xo, yo, zo, imgr, imgo, w, h, l, n0, n1, n2);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel linear interpolation 3d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel linear interpolation 3d" );
};

template <typename type>
void cuda_kernel_gradientx( std::vector<int> & grid, std::vector<int> & block,
                            const type * imgr, type * imgo, 
                            int n0, int n1, int n2)
{
    dim3 grd;
    dim3 blk;
    if (block[2] == 0)
    {
        grd = dim3(grid[0],grid[1]);
        blk = dim3(block[0],block[1]);
    }
    else
    {
        grd = dim3(grid[0],grid[1],grid[2]);
        blk = dim3(block[0],block[1],block[2]);
    };

    kernel_gradientx<<<grd,blk>>>(imgr, imgo, n0, n1, n2);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel gradient x" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel gradient x" );
};

template <typename type>
void cuda_kernel_gradienty( std::vector<int> & grid, std::vector<int> & block,
                            const type * imgr, type * imgo, 
                            int n0, int n1, int n2)
{
    dim3 grd;
    dim3 blk;
    if (block[2] == 0)
    {
        grd = dim3(grid[0],grid[1]);
        blk = dim3(block[0],block[1]);
    }
    else
    {
        grd = dim3(grid[0],grid[1],grid[2]);
        blk = dim3(block[0],block[1],block[2]);
    };
    kernel_gradienty<<<grd,blk>>>(imgr, imgo, n0, n1, n2);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel gradient y" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel gradient y" );
};

template <typename type>
void cuda_kernel_gradientz( std::vector<int> & grid, std::vector<int> & block,
                            const type * imgr, type * imgo, 
                            int n0, int n1, int n2)
{
    dim3 grd(grid[0],grid[1],grid[2]);
    dim3 blk(block[0],block[1],block[2]);
    kernel_gradientz<<<grd,blk>>>(imgr, imgo, n0, n1, n2);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel gradient z" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel gradient z" );
};

template <typename type>
void cuda_kernel_convolution_2d( std::vector<int> & grid, std::vector<int> & block,
                                 const type * imgr, const type * kern, //kernel width
                                 type * imgo, int kwidth, int n0, int n1)
{
    dim3 grd(grid[0],grid[1]);
    dim3 blk(block[0],block[1]);
    kernel_convolution_2d<<<grd,blk>>>(imgr, kern, imgo, kwidth, n0, n1);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel convolution 2d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel convolution 2d" );
};

template <typename type>
void cuda_kernel_convolution_3d( std::vector<int> & grid, std::vector<int> & block,
                                 const type * imgr, const type * kern, //kernel width
                                 type * imgo, int kwidth, int n0, int n1, int n2)
{
    dim3 grd(grid[0],grid[1],grid[2]);
    dim3 blk(block[0],block[1],block[2]);
    kernel_convolution_3d<<<grd,blk>>>(imgr, kern, imgo, kwidth, n0, n1, n2);
    imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel convolution 3d" );
    imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel convolution 3d" );
};


// template <typename type>
// void cuda_kernel_( std::vector<int> & grid, std::vector<int> & block, 
//                      )
// {
//     dim3 grd(grid[0],grid[1],grid[2]);
//     dim3 blk(block[0],block[1],block[2]);
//     kernel_<<<grd,blk>>>();
//     imart_assert_kernel( cudaPeekAtLastError(), "Fail to run kernel" );
//     imart_assert_kernel( cudaDeviceSynchronize(), "Fail to sync kernel" );
// };



// ===========================================
// Explicit instanciation
// ===========================================

// CASTINGS
template void cuda_kernel_cast<float,double>( std::vector<int> & grid, std::vector<int> & block, 
                       const float * vin, double * vout, int n );

template void cuda_kernel_cast<double,float>( std::vector<int> & grid, std::vector<int> & block, 
                       const double * vin, float * vout, int n );


template void cuda_kernel_cast<int,float>( std::vector<int> & grid, std::vector<int> & block, 
                       const int * vin, float * vout, int n );

template void cuda_kernel_cast<float,int>( std::vector<int> & grid, std::vector<int> & block, 
                       const float * vin, int * vout, int n );

template void cuda_kernel_cast<int,double>( std::vector<int> & grid, std::vector<int> & block, 
                       const int * vin, double * vout, int n );

template void cuda_kernel_cast<double,int>( std::vector<int> & grid, std::vector<int> & block, 
                       const double * vin, int * vout, int n );


template void cuda_kernel_cast<float,unsigned short>( std::vector<int> & grid, std::vector<int> & block, 
                       const float * vin, unsigned short * vout, int n );

template void cuda_kernel_cast<unsigned short,float>( std::vector<int> & grid, std::vector<int> & block, 
                       const unsigned short * vin, float * vout, int n );

template void cuda_kernel_cast<double,unsigned short>( std::vector<int> & grid, std::vector<int> & block, 
                       const double * vin, unsigned short * vout, int n );

template void cuda_kernel_cast<unsigned short,double>( std::vector<int> & grid, std::vector<int> & block, 
                       const unsigned short * vin, double * vout, int n );


template void cuda_kernel_cast<float,unsigned char>( std::vector<int> & grid, std::vector<int> & block, 
                       const float * vin, unsigned char * vout, int n );

template void cuda_kernel_cast<unsigned char,float>( std::vector<int> & grid, std::vector<int> & block, 
                       const unsigned char * vin, float * vout, int n );

template void cuda_kernel_cast<double,unsigned char>( std::vector<int> & grid, std::vector<int> & block, 
                       const double * vin, unsigned char * vout, int n );

template void cuda_kernel_cast<unsigned char,double>( std::vector<int> & grid, std::vector<int> & block, 
                       const unsigned char * vin, double * vout, int n );




template void cuda_kernel_assign<float>( std::vector<int> & grid, std::vector<int> & block,
                                         float * vin, float value, int n );

template void cuda_kernel_copy<float>( std::vector<int> & grid, std::vector<int> & block,
                                       const float * vin, float * vout, int n );

template void cuda_kernel_add<float>( std::vector<int> & grid, std::vector<int> & block,
                                      const float * vin1, const float * vin2, float * vout, int n );

template void cuda_kernel_sub<float>( std::vector<int> & grid, std::vector<int> & block,
                                      const float * vin1, const float * vin2, float * vout, int n );

template void cuda_kernel_mul<float>( std::vector<int> & grid, std::vector<int> & block,
                                      const float * vin1, const float * vin2, float * vout, int n );

template void cuda_kernel_div<float>( std::vector<int> & grid, std::vector<int> & block,
                                      const float * vin1, const float * vin2, float * vout, int n );

template void cuda_kernel_pow<float>( std::vector<int> & grid, std::vector<int> & block,
                                      const float * vin1, const float * vin2, float * vout, int n );

template void cuda_kernel_add_scalar<float>( std::vector<int> & grid, std::vector<int> & block, 
                                             const float * vin, float * vout, float scalar, int n );

template void cuda_kernel_sub_scalar<float>( std::vector<int> & grid, std::vector<int> & block, 
                                             const float * vin, float * vout, float scalar, int n );

template void cuda_kernel_sub_scalar_inv<float>( std::vector<int> & grid, std::vector<int> & block, 
                                             const float * vin, float * vout, float scalar, int n );

template void cuda_kernel_mul_scalar<float>( std::vector<int> & grid, std::vector<int> & block, 
                                             const float * vin, float * vout, float scalar, int n );

template void cuda_kernel_div_scalar<float>( std::vector<int> & grid, std::vector<int> & block, 
                                             const float * vin, float * vout, float scalar, int n );

template void cuda_kernel_div_scalar_inv<float>( std::vector<int> & grid, std::vector<int> & block, 
                                             const float * vin, float * vout, float scalar, int n );

template void cuda_kernel_pow_scalar<float>( std::vector<int> & grid, std::vector<int> & block, 
                                             const float * vin, float * vout, float scalar, int n );

template void cuda_kernel_pow_scalar_inv<float>( std::vector<int> & grid, std::vector<int> & block, 
                                             const float * vin, float * vout, float scalar, int n );

template void cuda_kernel_sum<float>( std::vector<int> & grid, std::vector<int> & block, const float * vin, float * vout, int n );

template void cuda_kernel_min<float>( std::vector<int> & grid, std::vector<int> & block, const float * vin, float * vout, int n );

template void cuda_kernel_max<float>( std::vector<int> & grid, std::vector<int> & block, const float * vin, float * vout, int n );


template void cuda_kernel_pad_2d<float>( std::vector<int> & grid, std::vector<int> & block, 
                            const float * vin, float * vout, int start0, int start1, 
                            int end0, int end1, int n0, int n1 );

template void cuda_kernel_unpad_2d<float>( std::vector<int> & grid, std::vector<int> & block, 
                            const float * vin, float * vout, int start0, int start1,
                            int end0, int end1, int n0, int n1);

template void cuda_kernel_pad_3d<float>( std::vector<int> & grid, std::vector<int> & block, 
                         const float * vin, float * vout, int start0, int start1, int start2,
                         int end0, int end1, int end2, int n0, int n1, int n2);

template void cuda_kernel_unpad_3d<float>( std::vector<int> & grid, std::vector<int> & block, 
                            const float * vin, float * vout, int start0, int start1, int start2,
                            int end0, int end1, int end2, int n0, int n1, int n2);

template void cuda_kernel_grid_2d<float>( std::vector<int> & grid, std::vector<int> & block,
                          float * x, float * y, double * sod, 
                          int n0, int n1);

template void cuda_kernel_grid_3d<float>( std::vector<int> & grid, std::vector<int> & block,
                          float * x, float * y, float * z, double * sod, 
                          int n0, int n1, int n2);

template void cuda_kernel_affine_2d<float>( std::vector<int> & grid, std::vector<int> & block,
                            const float * xin, const float * yin, 
                            float * xout, float * yout, 
                            const float * param, int n);

template void cuda_kernel_affine_3d<float>( std::vector<int> & grid, std::vector<int> & block,
                            const float * xin, const float * yin, const float * zin,
                            float * xout, float * yout, float * zout,
                            const float * param, int n) ;

template void cuda_kernel_affine_sod_2d<float>( std::vector<int> & grid, std::vector<int> & block,
                                const float * xin, const float * yin,
                                float * xout, float * yout,
                                const double * sod, int n);

template void cuda_kernel_affine_sod_3d<float>( std::vector<int> & grid, std::vector<int> & block,
                                const float * xin, const float * yin, const float * zin,
                                float * xout, float * yout, float * zout,
                                const double * sod, int n);

template void cuda_kernel_dfield_2d<float>( std::vector<int> & grid, std::vector<int> & block,
                            const float * xin, const float * yin,   // grid coordinates
                            const float * x, const float * y,       // vector field
                            float * xout, float * yout, int n);

template void cuda_kernel_dfield_3d<float>( std::vector<int> & grid, std::vector<int> & block,
                            const float * xin, const float * yin, const float * zin, // grid coordinates
                            const float * x, const float * y, const float * z,       // vector field
                            float * xout, float * yout, float * zout, int n);

template void cuda_kernel_nearest_interpolation_2d<float>( std::vector<int> & grid, std::vector<int> & block,
                                           const float * xo, const float * yo,
                                           const float * imgr, float * imgo,
                                           int w, int h,   //img ref width and height
                                           int n0, int n1); //img out dims

template void cuda_kernel_nearest_interpolation_3d<float>( std::vector<int> & grid, std::vector<int> & block, 
                                           const float * xo, const float * yo, const float * zo, 
                                           const float * imgr, float * imgo,
                                           int w, int h, int l,    //img ref width, height and length
                                           int n0, int n1, int n2 );

template void cuda_kernel_linear_interpolation_2d<float>( std::vector<int> & grid, std::vector<int> & block,
                                          const float * xo, const float * yo,
                                          const float * imgr, float * imgo,
                                          int w, int h,   //img ref width and height
                                          int n0, int n1); //img out dims

template void cuda_kernel_linear_interpolation_3d<float>( std::vector<int> & grid, std::vector<int> & block,
                                          const float * xo, const float * yo, const float * zo,
                                          const float * imgr, float * imgo,
                                          int w, int h, int l, //img ref width, height and length
                                          int n0, int n1, int n2);

template void cuda_kernel_gradientx<float>( std::vector<int> & grid, std::vector<int> & block,
                            const float * imgr, float * imgo, 
                            int n0, int n1, int n2);

template void cuda_kernel_gradienty<float>( std::vector<int> & grid, std::vector<int> & block,
                            const float * imgr, float * imgo, 
                            int n0, int n1, int n2);

template void cuda_kernel_gradientz<float>( std::vector<int> & grid, std::vector<int> & block,
                            const float * imgr, float * imgo, 
                            int n0, int n1, int n2);

template void cuda_kernel_convolution_2d<float>( std::vector<int> & grid, std::vector<int> & block,
                                 const float * imgr, const float * kern, //kernel width
                                 float * imgo, int kwidth, int n0, int n1);

template void cuda_kernel_convolution_3d<float>( std::vector<int> & grid, std::vector<int> & block,
                                 const float * imgr, const float * kern, //kernel width
                                 float * imgo, int kwidth, int n0, int n1, int n2);






template void cuda_kernel_assign<double>( std::vector<int> & grid, std::vector<int> & block,
                                         double * vin, double value, int n );

template void cuda_kernel_copy<double>( std::vector<int> & grid, std::vector<int> & block,
                                       const double * vin, double * vout, int n );

template void cuda_kernel_add<double>( std::vector<int> & grid, std::vector<int> & block,
                                      const double * vin1, const double * vin2, double * vout, int n );

template void cuda_kernel_sub<double>( std::vector<int> & grid, std::vector<int> & block,
                                      const double * vin1, const double * vin2, double * vout, int n );

template void cuda_kernel_mul<double>( std::vector<int> & grid, std::vector<int> & block,
                                      const double * vin1, const double * vin2, double * vout, int n );

template void cuda_kernel_div<double>( std::vector<int> & grid, std::vector<int> & block,
                                      const double * vin1, const double * vin2, double * vout, int n );

template void cuda_kernel_pow<double>( std::vector<int> & grid, std::vector<int> & block,
                                      const double * vin1, const double * vin2, double * vout, int n );

template void cuda_kernel_add_scalar<double>( std::vector<int> & grid, std::vector<int> & block, 
                                             const double * vin, double * vout, double scalar, int n );

template void cuda_kernel_sub_scalar<double>( std::vector<int> & grid, std::vector<int> & block, 
                                             const double * vin, double * vout, double scalar, int n );

template void cuda_kernel_sub_scalar_inv<double>( std::vector<int> & grid, std::vector<int> & block, 
                                             const double * vin, double * vout, double scalar, int n );

template void cuda_kernel_mul_scalar<double>( std::vector<int> & grid, std::vector<int> & block, 
                                             const double * vin, double * vout, double scalar, int n );

template void cuda_kernel_div_scalar<double>( std::vector<int> & grid, std::vector<int> & block, 
                                             const double * vin, double * vout, double scalar, int n );

template void cuda_kernel_div_scalar_inv<double>( std::vector<int> & grid, std::vector<int> & block, 
                                             const double * vin, double * vout, double scalar, int n );

template void cuda_kernel_pow_scalar<double>( std::vector<int> & grid, std::vector<int> & block, 
                                             const double * vin, double * vout, double scalar, int n );

template void cuda_kernel_pow_scalar_inv<double>( std::vector<int> & grid, std::vector<int> & block, 
                                             const double * vin, double * vout, double scalar, int n );

template void cuda_kernel_sum<double>( std::vector<int> & grid, std::vector<int> & block, const double * vin, double * vout, int n );

template void cuda_kernel_min<double>( std::vector<int> & grid, std::vector<int> & block, const double * vin, double * vout, int n );

template void cuda_kernel_max<double>( std::vector<int> & grid, std::vector<int> & block, const double * vin, double * vout, int n );

template void cuda_kernel_pad_2d<double>( std::vector<int> & grid, std::vector<int> & block, 
                            const double * vin, double * vout, int start0, int start1, 
                            int end0, int end1, int n0, int n1 );

template void cuda_kernel_unpad_2d<double>( std::vector<int> & grid, std::vector<int> & block, 
                            const double * vin, double * vout, int start0, int start1,
                            int end0, int end1, int n0, int n1);

template void cuda_kernel_pad_3d<double>( std::vector<int> & grid, std::vector<int> & block, 
                         const double * vin, double * vout, int start0, int start1, int start2,
                         int end0, int end1, int end2, int n0, int n1, int n2);

template void cuda_kernel_unpad_3d<double>( std::vector<int> & grid, std::vector<int> & block, 
                            const double * vin, double * vout, int start0, int start1, int start2,
                            int end0, int end1, int end2, int n0, int n1, int n2);

template void cuda_kernel_grid_2d<double>( std::vector<int> & grid, std::vector<int> & block,
                          double * x, double * y, double * sod, 
                          int n0, int n1);

template void cuda_kernel_grid_3d<double>( std::vector<int> & grid, std::vector<int> & block,
                          double * x, double * y, double * z, double * sod, 
                          int n0, int n1, int n2);

template void cuda_kernel_affine_2d<double>( std::vector<int> & grid, std::vector<int> & block,
                            const double * xin, const double * yin, 
                            double * xout, double * yout, 
                            const double * param, int n);

template void cuda_kernel_affine_3d<double>( std::vector<int> & grid, std::vector<int> & block,
                            const double * xin, const double * yin, const double * zin,
                            double * xout, double * yout, double * zout,
                            const double * param, int n) ;

template void cuda_kernel_affine_sod_2d<double>( std::vector<int> & grid, std::vector<int> & block,
                                const double * xin, const double * yin,
                                double * xout, double * yout,
                                const double * sod, int n);

template void cuda_kernel_affine_sod_3d<double>( std::vector<int> & grid, std::vector<int> & block,
                                const double * xin, const double * yin, const double * zin,
                                double * xout, double * yout, double * zout,
                                const double * sod, int n);

template void cuda_kernel_dfield_2d<double>( std::vector<int> & grid, std::vector<int> & block,
                            const double * xin, const double * yin,   // grid coordinates
                            const double * x, const double * y,       // vector field
                            double * xout, double * yout, int n);

template void cuda_kernel_dfield_3d<double>( std::vector<int> & grid, std::vector<int> & block,
                            const double * xin, const double * yin, const double * zin, // grid coordinates
                            const double * x, const double * y, const double * z,       // vector field
                            double * xout, double * yout, double * zout, int n);

template void cuda_kernel_nearest_interpolation_2d<double>( std::vector<int> & grid, std::vector<int> & block,
                                           const double * xo, const double * yo,
                                           const double * imgr, double * imgo,
                                           int w, int h,   //img ref width and height
                                           int n0, int n1); //img out dims

template void cuda_kernel_nearest_interpolation_3d<double>( std::vector<int> & grid, std::vector<int> & block, 
                                           const double * xo, const double * yo, const double * zo, 
                                           const double * imgr, double * imgo,
                                           int w, int h, int l,    //img ref width, height and length
                                           int n0, int n1, int n2);

template void cuda_kernel_linear_interpolation_2d<double>( std::vector<int> & grid, std::vector<int> & block,
                                          const double * xo, const double * yo,
                                          const double * imgr, double * imgo,
                                          int w, int h,   //img ref width and height
                                          int n0, int n1); //img out dims

template void cuda_kernel_linear_interpolation_3d<double>( std::vector<int> & grid, std::vector<int> & block,
                                          const double * xo, const double * yo, const double * zo,
                                          const double * imgr, double * imgo,
                                          int w, int h, int l, //img ref width, height and length
                                          int n0, int n1, int n2);

template void cuda_kernel_gradientx<double>( std::vector<int> & grid, std::vector<int> & block,
                            const double * imgr, double * imgo, 
                            int n0, int n1, int n2);

template void cuda_kernel_gradienty<double>( std::vector<int> & grid, std::vector<int> & block,
                            const double * imgr, double * imgo, 
                            int n0, int n1, int n2);

template void cuda_kernel_gradientz<double>( std::vector<int> & grid, std::vector<int> & block,
                            const double * imgr, double * imgo, 
                            int n0, int n1, int n2);

template void cuda_kernel_convolution_2d<double>( std::vector<int> & grid, std::vector<int> & block,
                                 const double * imgr, const double * kern, //kernel width
                                 double * imgo, int kwidth, int n0, int n1);

template void cuda_kernel_convolution_3d<double>( std::vector<int> & grid, std::vector<int> & block,
                                 const double * imgr, const double * kern, //kernel width
                                 double * imgo, int kwidth, int n0, int n1, int n2);





template void cuda_kernel_assign<int>( std::vector<int> & grid, std::vector<int> & block,
                                         int * vin, int value, int n );

template void cuda_kernel_copy<int>( std::vector<int> & grid, std::vector<int> & block,
                                       const int * vin, int * vout, int n );

template void cuda_kernel_add<int>( std::vector<int> & grid, std::vector<int> & block,
                                      const int * vin1, const int * vin2, int * vout, int n );

template void cuda_kernel_sub<int>( std::vector<int> & grid, std::vector<int> & block,
                                      const int * vin1, const int * vin2, int * vout, int n );

template void cuda_kernel_mul<int>( std::vector<int> & grid, std::vector<int> & block,
                                      const int * vin1, const int * vin2, int * vout, int n );

template void cuda_kernel_div<int>( std::vector<int> & grid, std::vector<int> & block,
                                      const int * vin1, const int * vin2, int * vout, int n );

template void cuda_kernel_pow<int>( std::vector<int> & grid, std::vector<int> & block,
                                      const int * vin1, const int * vin2, int * vout, int n );

template void cuda_kernel_add_scalar<int>( std::vector<int> & grid, std::vector<int> & block, 
                                             const int * vin, int * vout, int scalar, int n );

template void cuda_kernel_sub_scalar<int>( std::vector<int> & grid, std::vector<int> & block, 
                                             const int * vin, int * vout, int scalar, int n );

template void cuda_kernel_sub_scalar_inv<int>( std::vector<int> & grid, std::vector<int> & block, 
                                             const int * vin, int * vout, int scalar, int n );

template void cuda_kernel_mul_scalar<int>( std::vector<int> & grid, std::vector<int> & block, 
                                             const int * vin, int * vout, int scalar, int n );

template void cuda_kernel_div_scalar<int>( std::vector<int> & grid, std::vector<int> & block, 
                                             const int * vin, int * vout, int scalar, int n );

template void cuda_kernel_div_scalar_inv<int>( std::vector<int> & grid, std::vector<int> & block, 
                                             const int * vin, int * vout, int scalar, int n );

template void cuda_kernel_pow_scalar<int>( std::vector<int> & grid, std::vector<int> & block, 
                                             const int * vin, int * vout, int scalar, int n );

template void cuda_kernel_pow_scalar_inv<int>( std::vector<int> & grid, std::vector<int> & block, 
                                             const int * vin, int * vout, int scalar, int n );

template void cuda_kernel_sum<int>( std::vector<int> & grid, std::vector<int> & block, const int * vin, int * vout, int n );

template void cuda_kernel_min<int>( std::vector<int> & grid, std::vector<int> & block, const int * vin, int * vout, int n );

template void cuda_kernel_max<int>( std::vector<int> & grid, std::vector<int> & block, const int * vin, int * vout, int n );


template void cuda_kernel_pad_2d<int>( std::vector<int> & grid, std::vector<int> & block, 
                            const int * vin, int * vout, int start0, int start1, 
                            int end0, int end1, int n0, int n1 );

template void cuda_kernel_unpad_2d<int>( std::vector<int> & grid, std::vector<int> & block, 
                            const int * vin, int * vout, int start0, int start1,
                            int end0, int end1, int n0, int n1 );

template void cuda_kernel_pad_3d<int>( std::vector<int> & grid, std::vector<int> & block, 
                         const int * vin, int * vout, int start0, int start1, int start2,
                         int end0, int end1, int end2, int n0, int n1, int n2 );

template void cuda_kernel_unpad_3d<int>( std::vector<int> & grid, std::vector<int> & block, 
                            const int * vin, int * vout, int start0, int start1, int start2,
                            int end0, int end1, int end2, int n0, int n1, int n2 );

template void cuda_kernel_grid_2d<int>( std::vector<int> & grid, std::vector<int> & block,
                          int * x, int * y, double * sod, 
                          int n0, int n1 );

template void cuda_kernel_grid_3d<int>( std::vector<int> & grid, std::vector<int> & block,
                          int * x, int * y, int * z, double * sod, 
                          int n0, int n1, int n2 );

template void cuda_kernel_affine_2d<int>( std::vector<int> & grid, std::vector<int> & block,
                            const int * xin, const int * yin, 
                            int * xout, int * yout, 
                            const int * param, int n );

template void cuda_kernel_affine_3d<int>( std::vector<int> & grid, std::vector<int> & block,
                            const int * xin, const int * yin, const int * zin,
                            int * xout, int * yout, int * zout,
                            const int * param, int n );

template void cuda_kernel_affine_sod_2d<int>( std::vector<int> & grid, std::vector<int> & block,
                                const int * xin, const int * yin,
                                int * xout, int * yout,
                                const double * sod, int n);

template void cuda_kernel_affine_sod_3d<int>( std::vector<int> & grid, std::vector<int> & block,
                                const int * xin, const int * yin, const int * zin,
                                int * xout, int * yout, int * zout,
                                const double * sod, int n );

template void cuda_kernel_dfield_2d<int>( std::vector<int> & grid, std::vector<int> & block,
                            const int * xin, const int * yin,   // grid coordinates
                            const int * x, const int * y,       // vector field
                            int * xout, int * yout, int n );

template void cuda_kernel_dfield_3d<int>( std::vector<int> & grid, std::vector<int> & block,
                            const int * xin, const int * yin, const int * zin, // grid coordinates
                            const int * x, const int * y, const int * z,       // vector field
                            int * xout, int * yout, int * zout, int n );

// template void cuda_kernel_nearest_interpolation_2d<int>( std::vector<int> & grid, std::vector<int> & block,
//                                            const int * xo, const int * yo,
//                                            const int * imgr, int * imgo,
//                                            int w, int h,   //img ref width and height
//                                            int n0, int n1); //img out dims

// template void cuda_kernel_nearest_interpolation_3d<int>( std::vector<int> & grid, std::vector<int> & block, 
//                                            const int * xo, const int * yo, const int * zo, 
//                                            const int * imgr, int * imgo,
//                                            int w, int h, int l,    //img ref width, height and length
//                                            int n0, int n1, int n2);

// template void cuda_kernel_linear_interpolation_2d<int>( std::vector<int> & grid, std::vector<int> & block,
//                                           const int * xo, const int * yo,
//                                           const int * imgr, int * imgo,
//                                           int w, int h,   //img ref width and height
//                                           int n0, int n1); //img out dims

// template void cuda_kernel_linear_interpolation_3d<int>( std::vector<int> & grid, std::vector<int> & block,
//                                           const int * xo, const int * yo, const int * zo,
//                                           const int * imgr, int * imgo,
//                                           int w, int h, int l, //img ref width, height and length
//                                           int n0, int n1, int n2);

template void cuda_kernel_gradientx<int>( std::vector<int> & grid, std::vector<int> & block,
                            const int * imgr, int * imgo, 
                            int n0, int n1, int n2);

template void cuda_kernel_gradienty<int>( std::vector<int> & grid, std::vector<int> & block,
                            const int * imgr, int * imgo, 
                            int n0, int n1, int n2);

template void cuda_kernel_gradientz<int>( std::vector<int> & grid, std::vector<int> & block,
                            const int * imgr, int * imgo, 
                            int n0, int n1, int n2);

template void cuda_kernel_convolution_2d<int>( std::vector<int> & grid, std::vector<int> & block,
                                 const int * imgr, const int * kern, //kernel width
                                 int * imgo, int kwidth, int n0, int n1);

template void cuda_kernel_convolution_3d<int>( std::vector<int> & grid, std::vector<int> & block,
                                 const int * imgr, const int * kern, //kernel width
                                 int * imgo, int kwidth, int n0, int n1, int n2);





template void cuda_kernel_assign<unsigned short>( std::vector<int> & grid, std::vector<int> & block,
                                         unsigned short * vin, unsigned short value, int n );

template void cuda_kernel_copy<unsigned short>( std::vector<int> & grid, std::vector<int> & block,
                                       const unsigned short * vin, unsigned short * vout, int n );


template void cuda_kernel_assign<unsigned int>( std::vector<int> & grid, std::vector<int> & block,
                                         unsigned int * vin, unsigned int value, int n );

template void cuda_kernel_copy<unsigned int>( std::vector<int> & grid, std::vector<int> & block,
                                       const unsigned int * vin, unsigned int * vout, int n );


template void cuda_kernel_assign<unsigned char>( std::vector<int> & grid, std::vector<int> & block,
                                         unsigned char * vin, unsigned char value, int n );

template void cuda_kernel_copy<unsigned char>( std::vector<int> & grid, std::vector<int> & block,
                                       const unsigned char * vin, unsigned char * vout, int n );


template void cuda_kernel_assign<short>( std::vector<int> & grid, std::vector<int> & block,
                                         short * vin, short value, int n );

template void cuda_kernel_copy<short>( std::vector<int> & grid, std::vector<int> & block,
                                       const short * vin, short * vout, int n );


template void cuda_kernel_assign<char>( std::vector<int> & grid, std::vector<int> & block,
                                         char * vin, char value, int n );

template void cuda_kernel_copy<char>( std::vector<int> & grid, std::vector<int> & block,
                                       const char * vin, char * vout, int n );
