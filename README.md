# lid-driven-cavity-2D

<p>Industrial CFD simulator for lid-driven cavity flow, implemented in Fortran with a modular architecture (params, grid, numerics, solver, I/O). Supports 2D uniform Cartesian meshes and second-order differential operators (gradient, divergence, Laplacian, curl) with advective schemes upwind/QUICK. Dirichlet/Neumann boundary conditions and adaptive CFL control ensure incompressibility and momentum conservation. The solver applies Chorin’s fractional step projection method, solving the Poisson equation via an optimized Red-Black SOR, with adaptive residual computation and OpenMP parallelization, maximizing cache efficiency and avoiding race conditions. Hierarchical HDF5 I/O with DEFLATE compression, optimized chunking, and full metadata guarantees bit-exact reproducibility, comprehensive logging, and traceability of physical and numerical parameters. Automatic diagnostics for kinetic energy, enstrophy, vorticity, and divergence, validation against Ghia et al. (1982) benchmarks, and real-time NaN/Inf detection ensure numerical robustness. Extensible and configurable design for RK2/RK4 schemes, multigrid/Krylov solvers, MPI, and AMR, optimizing scalability, throughput, and computational robustness for research, advanced education, and industrial development.</p>

<p>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-7D9EC0?logo=github" />
  </a>
  <img src="https://img.shields.io/badge/Docker-Active-29ABE2?logo=docker" />
  <img src="https://img.shields.io/badge/Fortran-2023-F05032?logo=gnu" />
  <img src="https://img.shields.io/badge/GNU_Fortran-14.2.0-E95420?logo=gnu" />
  <img src="https://img.shields.io/badge/CMake-3.28.3-064F8C?logo=cmake" />
  <img src="https://img.shields.io/badge/HDF5-1.14.3-4393D3?logo=hdf5" />
</p>

## 📂 Project Structure
```text
lid-driven-cavity-2D/
├── CMakeLists.txt          # Build configuration
├── docker-compose.yml      # Docker compose setup
├── Dockerfile              # Docker container definition
├── LICENSE                 # Project license
├── README.md               # Project documentation
├── build/                  # Build output directory (generated)
└── src/
    └── fortran/
        ├── grid_mod.f90       # Mesh generation and grid utilities
        ├── io_hdf5_mod.f90    # HDF5 I/O routines (read/write, metadata, compression)
        ├── main.f90           # Main program entry point
        ├── numerics_mod.f90   # Differential operators, advective schemes, validation
        ├── params.f90.in      # Template for project parameters (filled by CMake)
        └── solver_mod.f90     # Fractional step solver, SOR Red-Black, diagnostics
```
