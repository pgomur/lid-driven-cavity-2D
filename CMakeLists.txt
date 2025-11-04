cmake_minimum_required(VERSION 3.26)
project(lid_driven_cavity_2D Fortran)

set(CMAKE_Fortran_STANDARD 2023)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_FLAGS "-O2 -Wall -fopenmp")

# ====================================
# Physical and numerical parameters (overwriteable)
# ====================================
set(PARAM_NX 64 CACHE STRING "Número de celdas en X")
set(PARAM_NY 64 CACHE STRING "Número de celdas en Y")
set(PARAM_LX 1.0 CACHE STRING "Largo del dominio X")
set(PARAM_LY 1.0 CACHE STRING "Largo del dominio Y")
set(PARAM_XMIN 0.0 CACHE STRING "X mínimo")
set(PARAM_XMAX 1.0 CACHE STRING "X máximo")
set(PARAM_YMIN 0.0 CACHE STRING "Y mínimo")
set(PARAM_YMAX 1.0 CACHE STRING "Y máximo")

set(PARAM_RE 1000.0 CACHE STRING "Reynolds")
set(PARAM_ULID 1.0 CACHE STRING "Velocidad del lid")
set(PARAM_RHO 1.0 CACHE STRING "Densidad")
set(PARAM_NU 0.001 CACHE STRING "Viscosidad")
set(PARAM_GRAVITY 0.0 CACHE STRING "Gravedad")
set(PARAM_BETA 0.0 CACHE STRING "Beta Boussinesq")

set(PARAM_T_END 1.0 CACHE STRING "Tiempo final")
set(PARAM_DT 0.0024 CACHE STRING "Paso de tiempo")
set(PARAM_CFL 0.5 CACHE STRING "CFL")
set(PARAM_TIME_SCHEME "RK4" CACHE STRING "Esquema temporal")

# Output / post-processing
set(PARAM_OUTPUT_EVERY 10 CACHE STRING "Frecuencia salida")
set(PARAM_OUTPUT_FORMAT "hdf5" CACHE STRING "Formato salida")
set(PARAM_OUTPUT_DIR "results" CACHE STRING "Directorio salida")
set(PARAM_SAVE_VELOCITY TRUE CACHE STRING "Guardar velocidad")
set(PARAM_SAVE_PRESSURE FALSE CACHE STRING "Guardar presión")

# Boundary conditions
set(PARAM_BC_TOP "Dirichlet" CACHE STRING "BC top")
set(PARAM_BC_BOTTOM "Dirichlet" CACHE STRING "BC bottom")
set(PARAM_BC_LEFT "Dirichlet" CACHE STRING "BC left")
set(PARAM_BC_RIGHT "Dirichlet" CACHE STRING "BC right")

# Solver
set(PARAM_SOLVER_TYPE "explicit" CACHE STRING "Tipo de solver")
set(PARAM_MAX_ITERS 10000 CACHE STRING "Número máximo de iteraciones")
set(PARAM_TOL 1.0e-8 CACHE STRING "Tolerancia solver")
set(PARAM_LINEAR_SOLVER "CG" CACHE STRING "Solver lineal")
set(PARAM_PRECONDITIONER "ILU" CACHE STRING "Precondicionador")

# Discretization and numerical
set(PARAM_NUM_SCHEME "central" CACHE STRING "Esquema numérico")
set(PARAM_ORDER 2 CACHE STRING "Orden del esquema")
set(PARAM_STABILIZATION 0.0 CACHE STRING "Estabilización")

# Debug
set(PARAM_VERBOSE TRUE CACHE STRING "Verbose")
set(PARAM_CHECK_BOUNDS TRUE CACHE STRING "Comprobar límites")

# ====================================
# Generate params.f90
# ====================================
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/fortran/params.f90.in
               ${CMAKE_BINARY_DIR}/generated_params.f90
               @ONLY)

# ====================================
# Modules and directories.
# ====================================
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
include_directories(${CMAKE_BINARY_DIR})

# Source files
set(SOURCES
  ${CMAKE_BINARY_DIR}/generated_params.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fortran/io_hdf5_mod.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fortran/solver_mod.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fortran/numerics_mod.f90 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fortran/grid_mod.f90 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fortran/main.f90
)

# Executable
add_executable(lid_driven_cavity ${SOURCES})
set_target_properties(lid_driven_cavity PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# HDF5 Fortran
set(HDF5_ROOT /opt/hdf5)
find_package(HDF5 REQUIRED COMPONENTS Fortran)
if(HDF5_Fortran_FOUND)
    message(STATUS "HDF5 Fortran found")
    target_include_directories(lid_driven_cavity PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
    target_link_libraries(lid_driven_cavity PRIVATE ${HDF5_Fortran_LIBRARIES})
else()
    message(FATAL_ERROR "HDF5 Fortran no encontrado")
endif()

# OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(lid_driven_cavity PRIVATE OpenMP::OpenMP_Fortran)
endif()
