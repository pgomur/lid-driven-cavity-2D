# ===============================
# Optimized Dockerfile for Fortran + HDF5
# ===============================

# Base image
FROM ubuntu:24.04

# ===============================
# Environment variables for HDF5 and GCC
# ===============================
ENV HDF5_VERSION=1.14.3
ENV HDF5_DIR=/opt/hdf5
ENV PATH=${HDF5_DIR}/bin:$PATH
ENV LD_LIBRARY_PATH=${HDF5_DIR}/lib:$LD_LIBRARY_PATH
ENV LIBRARY_PATH=${HDF5_DIR}/lib:$LIBRARY_PATH
ENV CPATH=${HDF5_DIR}/include:$CPATH

# ===============================
# Install basic dependencies
# ===============================
RUN apt-get update && apt-get install -y --no-install-recommends \
        software-properties-common \
        wget \
        build-essential \
        cmake \
        m4 \
        make \
        autoconf \
        automake \
        libtool \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ===============================
# Install GCC/GFortran 14
# ===============================
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y && apt-get update && \
    apt-get install -y --no-install-recommends \
        gcc-14 g++-14 gfortran-14 \
    && update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-14 100 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ===============================
# Install HDF5 with Fortran support
# ===============================
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz && \
    tar -xzf hdf5-${HDF5_VERSION}.tar.gz && \
    cd hdf5-${HDF5_VERSION} && \
    ./configure --prefix=${HDF5_DIR} --enable-fortran --enable-shared && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf hdf5-${HDF5_VERSION} hdf5-${HDF5_VERSION}.tar.gz

# ===============================
# Set working directory inside the container
# ===============================
WORKDIR /app

# ===============================
# Copy CMakeLists.txt and Fortran source code
# ===============================
COPY CMakeLists.txt .
COPY src/ ./src

# ===============================
# Prepare build folder and compile Fortran project
# ===============================
RUN rm -rf build && mkdir -p build && cd build && cmake .. && make -j$(nproc)

# ===============================
# Default command when the container starts
# ===============================
CMD ["bash"]
