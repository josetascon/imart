# IMAge RegisTration Library (IMART)

This is a registration library using cpu, and gpu.
The aim is to achieve real-time deformable registration and reduce computational time for deformable registration.

## Desing decisions:

* Library is header only. This simplify template creation. Template class definition short and only brief descriptions. Functions names intuitive.
* ITK core for image and transformations interface (reading and writing files).
* A container is defined as the main class to support vector operations. These classes are: vector_cpu, vector_opencl, vector_cuda.

## Dependencies:
* ITK (read/write images)
* VTK (visualize images)
* Boost (program options)
* OpenMP, OpenCL, CUDA

### Ubuntu

Installing all dependencies in ubuntu.\

General dependencies:\
\# apt install -y cmake git wget\
OpenCL dependencies:\
\# apt install -y opencl-headers opencl-c-headers opencl-clhpp-headers\
\# apt install -y clinfo ocl-icd-libopencl1 ocl-icd-opencl-dev beignet-opencl-icd\
Boost libraries:\
\# apt install -y libboost1.65-dev libboost-program-options1.65-dev\
ITK libraries:\
\# apt install -y libinsighttoolkit4.12 libinsighttoolkit4-dev insighttoolkit4-examples\
VTK libraries:\
\# apt install -y tzdata libvtk7.1 libvtk7-dev vtk7 vtk7-examples

### Arch Linux

Install all dependencies.\

\# pacman -S cmake git opencl-headers cuda boost vtk

Install form AUR packages:\
insight-toolkit

## Build

In the project folder use cmake to build.

mkdir build\
cd build\
cmake .. -DBUILD_EXAMPLES=ON -DBUILD_TESTS=ON -DBUILD_BENCHMARKS=ON -DIMART_WITH_OPENCL=ON -DIMART_WITH_CUDA=ON -DIMART_WITH_CLFFT=ON

Change ON or OFF depending on your dependencies.