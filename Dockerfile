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



# docker-compose build
# docker-compose up -d
# docker-compose up --build 

# cd build            # tu directorio de compilación
# cmake ..            # genera Makefiles
# cmake --build .     # compila
# ./bin/lid_driven_cavity  # ejecuta el binario final



# cd build && cmake .. -DHDF5_ROOT=/opt/hdf5 -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O2 -Wall -fopenmp" && cmake --build . -j$(nproc) && ./bin/lid_driven_cavity

# cd build && cmake .. -DHDF5_ROOT=/opt/hdf5 -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O2 -Wall -fopenmp" -DPARAM_NX=206 -DPARAM_NY=206 -DPARAM_LX=1.0 -DPARAM_LY=1.0 -DPARAM_XMIN=0.0 -DPARAM_XMAX=1.0 -DPARAM_YMIN=0.0 -DPARAM_YMAX=1.0 -DPARAM_RE=1000.0 -DPARAM_ULID=1.0 -DPARAM_RHO=1.0 -DPARAM_NU=0.001 -DPARAM_GRAVITY=0.0 -DPARAM_BETA=0.0 -DPARAM_T_END=1.0 -DPARAM_DT=0.0024 -DPARAM_CFL=0.5 -DPARAM_TIME_SCHEME=RK4 -DPARAM_OUTPUT_EVERY=10 -DPARAM_OUTPUT_FORMAT=hdf5 -DPARAM_OUTPUT_DIR=results -DPARAM_SAVE_VELOCITY=TRUE -DPARAM_SAVE_PRESSURE=FALSE -DPARAM_BC_TOP=Dirichlet -DPARAM_BC_BOTTOM=Dirichlet -DPARAM_BC_LEFT=Dirichlet -DPARAM_BC_RIGHT=Dirichlet -DPARAM_SOLVER_TYPE=explicit -DPARAM_MAX_ITERS=10000 -DPARAM_TOL=1.0e-8 -DPARAM_LINEAR_SOLVER=CG -DPARAM_PRECONDITIONER=ILU -DPARAM_NUM_SCHEME=central -DPARAM_ORDER=2 -DPARAM_STABILIZATION=0.0 -DPARAM_VERBOSE=TRUE -DPARAM_CHECK_BOUNDS=TRUE && cmake --build . -j$(nproc) && ./bin/lid_driven_cavity
# cd build && cmake .. -DHDF5_ROOT=/opt/hdf5 -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops -fopenmp -Wall" -DPARAM_NX=64 -DPARAM_NY=64 -DPARAM_LX=1.0 -DPARAM_LY=1.0 -DPARAM_XMIN=0.0 -DPARAM_XMAX=1.0 -DPARAM_YMIN=0.0 -DPARAM_YMAX=1.0 -DPARAM_RE=1000.0 -DPARAM_ULID=1.0 -DPARAM_RHO=1.0 -DPARAM_NU=0.001 -DPARAM_GRAVITY=0.0 -DPARAM_BETA=0.0 -DPARAM_T_END=1.0 -DPARAM_DT=0.0024 -DPARAM_CFL=0.5 -DPARAM_TIME_SCHEME=RK4 -DPARAM_OUTPUT_EVERY=10 -DPARAM_OUTPUT_FORMAT=hdf5 -DPARAM_OUTPUT_DIR=results -DPARAM_SAVE_VELOCITY=TRUE -DPARAM_SAVE_PRESSURE=FALSE -DPARAM_BC_TOP=Dirichlet -DPARAM_BC_BOTTOM=Dirichlet -DPARAM_BC_LEFT=Dirichlet -DPARAM_BC_RIGHT=Dirichlet -DPARAM_SOLVER_TYPE=explicit -DPARAM_MAX_ITERS=10000 -DPARAM_TOL=1.0e-8 -DPARAM_LINEAR_SOLVER=CG -DPARAM_PRECONDITIONER=ILU -DPARAM_NUM_SCHEME=central -DPARAM_ORDER=2 -DPARAM_STABILIZATION=0.0 -DPARAM_VERBOSE=TRUE -DPARAM_CHECK_BOUNDS=TRUE && cmake --build . -j$(nproc) && export OMP_NUM_THREADS=$(nproc) && ./bin/lid_driven_cavity

# cd build && cmake .. -DHDF5_ROOT=/opt/hdf5 -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops -fopenmp -Wall" -DPARAM_NX=32 -DPARAM_NY=32 -DPARAM_LX=1.0 -DPARAM_LY=1.0 -DPARAM_XMIN=0.0 -DPARAM_XMAX=1.0 -DPARAM_YMIN=0.0 -DPARAM_YMAX=1.0 -DPARAM_RE=100.0 -DPARAM_ULID=1.0 -DPARAM_RHO=1.0 -DPARAM_NU=0.01 -DPARAM_GRAVITY=0.0 -DPARAM_BETA=0.0 -DPARAM_T_END=1.0 -DPARAM_DT=0.005 -DPARAM_CFL=0.5 -DPARAM_TIME_SCHEME=RK4 -DPARAM_OUTPUT_EVERY=20 -DPARAM_OUTPUT_FORMAT=hdf5 -DPARAM_OUTPUT_DIR=results -DPARAM_SAVE_VELOCITY=TRUE -DPARAM_SAVE_PRESSURE=TRUE -DPARAM_BC_TOP=Dirichlet -DPARAM_BC_BOTTOM=Dirichlet -DPARAM_BC_LEFT=Dirichlet -DPARAM_BC_RIGHT=Dirichlet -DPARAM_SOLVER_TYPE=explicit -DPARAM_MAX_ITERS=5000 -DPARAM_TOL=1.0e-6 -DPARAM_LINEAR_SOLVER=CG -DPARAM_PRECONDITIONER=ILU -DPARAM_NUM_SCHEME=central -DPARAM_ORDER=2 -DPARAM_STABILIZATION=0.0 -DPARAM_VERBOSE=TRUE -DPARAM_CHECK_BOUNDS=TRUE && cmake --build . -j$(nproc) && export OMP_NUM_THREADS=$(nproc) && ./bin/lid_driven_cavity


# {
#         "PARAM_NX": 206,
#         "PARAM_NY": 206,
#         "PARAM_LX": 1.0,
#         "PARAM_LY": 1.0,
#         "PARAM_XMIN": 0.0,
#         "PARAM_XMAX": 1.0,
#         "PARAM_YMIN": 0.0,
#         "PARAM_YMAX": 1.0,
#         "PARAM_RE": 100.0,
#         "PARAM_ULID": 1.0,
#         "PARAM_RHO": 1.0,
#         "PARAM_NU": 0.01,
#         "PARAM_GRAVITY": 0.0,
#         "PARAM_BETA": 0.0,
#         "PARAM_T_END": 0.5,
#         "PARAM_DT": 0.002,
#         "PARAM_CFL": 0.5,
#         "PARAM_TIME_SCHEME": "RK4",
#         "PARAM_OUTPUT_EVERY": 20,
#         "PARAM_OUTPUT_FORMAT": "hdf5",
#         "PARAM_OUTPUT_DIR": "results",
#         "PARAM_SAVE_VELOCITY": "TRUE",
#         "PARAM_SAVE_PRESSURE": "FALSE",
#         "PARAM_BC_TOP": "Dirichlet",
#         "PARAM_BC_BOTTOM": "Dirichlet",
#         "PARAM_BC_LEFT": "Dirichlet",
#         "PARAM_BC_RIGHT": "Dirichlet",
#         "PARAM_SOLVER_TYPE": "explicit",
#         "PARAM_MAX_ITERS": 5000,
#         "PARAM_TOL": 1e-6,
#         "PARAM_LINEAR_SOLVER": "CG",
#         "PARAM_PRECONDITIONER": "ILU",
#         "PARAM_NUM_SCHEME": "central",
#         "PARAM_ORDER": 2,
#         "PARAM_STABILIZATION": 0.0,
#         "PARAM_VERBOSE": "TRUE",
#         "PARAM_CHECK_BOUNDS": "TRUE"
#      }