FROM nvidia/cuda
LABEL maintainer="josedtascon@gmail.com"

# =========================================
#       Install additional packages
# =========================================

RUN apt update

# Install build tools
RUN apt install -y cmake git wget

# Install opencl
RUN apt install -y opencl-headers opencl-c-headers opencl-clhpp-headers
RUN apt install -y clinfo ocl-icd-libopencl1 ocl-icd-opencl-dev beignet-opencl-icd

# Install boost
# RUN apt install libboost-all-dev
RUN apt install -y libboost1.65-dev libboost-program-options1.65-dev

# Install tzdata, required later for vtk
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Berlin
RUN apt install -y tzdata

RUN apt install -y libinsighttoolkit4.12 libinsighttoolkit4-dev insighttoolkit4-examples
RUN apt install -y libvtk7.1 libvtk7-dev vtk7 vtk7-examples

# =========================================
#               Configure
# =========================================
RUN mkdir -p /home/project/imart/

# Change later for git clone ...
COPY imart.tar.bz2 /home/project/imart/

#RUN cd /home/project/imart/
#RUN tar -xvf imart.tar.bz2
#RUN mkdir -p build
#RUN cd build
#RUN cmake ..