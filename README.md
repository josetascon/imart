IMAge RegisTration Library (IMART)

This is a registration library using cpu, and gpu.
The aim is to achieve real-time deformable registration and reduce computational time for deformable registration.

Dependencies:
* ITK (read/write images)
* VTK (visualize images)
* Boost (program options)
* OpenMP, OpenCL, CUDA

Desing decisions:

* Library is header only. This simplify template creation. Template class definition short and only brief descriptions. Functions names intuitive.
* ITK core for image and transformations interface (reading and writing files).
* A container is defined as the main class to support vector operations. These classes are: vector_cpu, vector_opencl, vector_cuda.



