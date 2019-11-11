FROM nvidia/cuda
LABEL maintainer="josedtascon@gmail.com"

# =========================================
#       Install additional packages
# =========================================

RUN apt-get update
RUN apt-get install -y \
    wget \
    git \
    cmake \ 
    clinfo \
    opencl-headers \
    opencl-c-headers \
    opencl-clhpp-headers \
    ocl-icd-libopencl1 \
    ocl-icd-opencl-dev \
    beignet-opencl-icd

# =========================================
#               Configure
# =========================================
RUN mkdir /home/project