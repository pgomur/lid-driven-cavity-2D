# lid-driven-cavity-2D

<p>Industrial CFD simulator for lid-driven cavity flow, implemented in Fortran with a modular architecture (params, grid, numerics, solver, I/O). Supports 2D uniform Cartesian meshes and second-order differential operators (gradient, divergence, Laplacian, curl) with advective schemes upwind/QUICK. Dirichlet/Neumann boundary conditions and adaptive CFL control ensure incompressibility and momentum conservation. The solver applies Chorin’s fractional step projection method, solving the Poisson equation via an optimized Red-Black SOR, with adaptive residual computation and OpenMP parallelization, maximizing cache efficiency and avoiding race conditions. Hierarchical HDF5 I/O with DEFLATE compression, optimized chunking, and full metadata guarantees bit-exact reproducibility, comprehensive logging, and traceability of physical and numerical parameters. Automatic diagnostics for kinetic energy, enstrophy, vorticity, and divergence, validation against Ghia et al. (1982) benchmarks, and real-time NaN/Inf detection ensure numerical robustness. Extensible and configurable design for RK2/RK4 schemes, multigrid/Krylov solvers, MPI, and AMR, optimizing scalability, throughput, and computational robustness for research, advanced education, and industrial development.</p>

<p>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-7D9EC0?logo=github" /></a>
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

<h3>🚀 Running the Project</h3>

<h4>Build the Docker image</h4>
<pre><code class="bash">docker-compose build</code></pre>

<h4>Start the container</h4>
<pre><code class="bash">docker-compose up -d</code></pre>

<h4>Enter the container shell</h4>
<pre><code class="bash">docker exec -it lid_driven_cavity_2d_container bash</code></pre>

<h4>Compile and run the simulator inside the container</h4>
<p>Inside the container, run the following command to compile and execute the CFD simulator:</p>
<pre><code class="bash">cd build && cmake .. -DHDF5_ROOT=/opt/hdf5 -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O3 -march=native -funroll-loops -fopenmp -Wall" -DPARAM_NX=32 -DPARAM_NY=32 -DPARAM_LX=1.0 -DPARAM_LY=1.0 -DPARAM_XMIN=0.0 -DPARAM_XMAX=1.0 -DPARAM_YMIN=0.0 -DPARAM_YMAX=1.0 -DPARAM_RE=100.0 -DPARAM_ULID=1.0 -DPARAM_RHO=1.0 -DPARAM_NU=0.01 -DPARAM_GRAVITY=0.0 -DPARAM_BETA=0.0 -DPARAM_T_END=1.0 -DPARAM_DT=0.005 -DPARAM_CFL=0.5 -DPARAM_TIME_SCHEME=RK4 -DPARAM_OUTPUT_EVERY=20 -DPARAM_OUTPUT_FORMAT=hdf5 -DPARAM_OUTPUT_DIR=results -DPARAM_SAVE_VELOCITY=TRUE -DPARAM_SAVE_PRESSURE=TRUE -DPARAM_BC_TOP=Dirichlet -DPARAM_BC_BOTTOM=Dirichlet -DPARAM_BC_LEFT=Dirichlet -DPARAM_BC_RIGHT=Dirichlet -DPARAM_SOLVER_TYPE=explicit -DPARAM_MAX_ITERS=5000 -DPARAM_TOL=1.0e-6 -DPARAM_LINEAR_SOLVER=CG -DPARAM_PRECONDITIONER=ILU -DPARAM_NUM_SCHEME=central -DPARAM_ORDER=2 -DPARAM_STABILIZATION=0.0 -DPARAM_VERBOSE=TRUE -DPARAM_CHECK_BOUNDS=TRUE && cmake --build . -j$(nproc) && export OMP_NUM_THREADS=$(nproc) && ./bin/lid_driven_cavity</code></pre>

<p>💡 <strong>Note:</strong> This command will build the project using all CPU cores and then run the simulator with the configured parameters.</p>
<p>📂 <strong>Simulation Results:</strong> The output files are saved inside the <code>build/results</code> folder in <strong>HDF5</strong> format.</p>
<h4>📊 HDF5 Dataset Structure</h4>
<pre><code>/lid_cavity.h5
├── x_coords [1D array]        (X coordinates)
├── y_coords [1D array]        (Y coordinates)
├── U_000000 [2D array]        (X velocity, step 0)
├── V_000000 [2D array]        (Y velocity, step 0)
├── P_000000 [2D array]        (Pressure, step 0)
├── time_000000 [attribute]    (Physical time)
├── step_000000 [attribute]    (Step number)
├── ...                        (Subsequent steps)
└── Attributes                 (Re, nu, dt, version, etc.)
</code></pre>

<h4>Explanation of the command:</h4>
<ul>
  <li><code>cd build</code><br>
      Move into the <strong>build directory</strong>, where the project will be compiled.</li>
  <li><code>cmake ..</code><br>
      Configure the project using <strong>CMake</strong>, pointing to the source code in the parent folder. This prepares the Makefiles or build system for compilation.</li>
  <li><code>cmake --build . -j$(nproc)</code><br>
      Compile the project using all available CPU cores (<code>-j$(nproc)</code>). Produces the executable in the <code>bin/</code> folder.</li>
  <li><code>export OMP_NUM_THREADS=$(nproc)</code><br>
      Set the number of threads for <strong>OpenMP parallelization</strong> to use all CPU cores.</li>
  <li><code>./bin/lid_driven_cavity</code><br>
      Run the <strong>simulator executable</strong>, which will start the CFD simulation using the configured parameters.</li>
</ul>

<h4>⚙️ Simulation Parameters</h4>
<table>
  <thead>
    <tr>
      <th>Parameter</th>
      <th>Type</th>
      <th>Typical Values</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr><td><code>PARAM_NX</code></td><td>int</td><td>&gt;1</td><td>Number of nodes in X-direction</td></tr>
    <tr><td><code>PARAM_NY</code></td><td>int</td><td>&gt;1</td><td>Number of nodes in Y-direction</td></tr>
    <tr><td><code>PARAM_LX</code></td><td>real</td><td>&gt;0</td><td>Domain length in X</td></tr>
    <tr><td><code>PARAM_LY</code></td><td>real</td><td>&gt;0</td><td>Domain length in Y</td></tr>
    <tr><td><code>PARAM_XMIN</code></td><td>real</td><td>&ge;0</td><td>Minimum X coordinate</td></tr>
    <tr><td><code>PARAM_XMAX</code></td><td>real</td><td>&gt; XMIN</td><td>Maximum X coordinate</td></tr>
    <tr><td><code>PARAM_YMIN</code></td><td>real</td><td>&ge;0</td><td>Minimum Y coordinate</td></tr>
    <tr><td><code>PARAM_YMAX</code></td><td>real</td><td>&gt; YMIN</td><td>Maximum Y coordinate</td></tr>
    <tr><td><code>PARAM_RE</code></td><td>real</td><td>&gt;0</td><td>Reynolds number</td></tr>
    <tr><td><code>PARAM_ULID</code></td><td>real</td><td>&ge;0</td><td>Lid velocity</td></tr>
    <tr><td><code>PARAM_RHO</code></td><td>real</td><td>&gt;0</td><td>Fluid density</td></tr>
    <tr><td><code>PARAM_NU</code></td><td>real</td><td>&gt;0</td><td>Kinematic viscosity ν</td></tr>
    <tr><td><code>PARAM_GRAVITY</code></td><td>real</td><td>&ge;0</td><td>Gravitational acceleration (Boussinesq)</td></tr>
    <tr><td><code>PARAM_BETA</code></td><td>real</td><td>&ge;0</td><td>Thermal expansion coefficient (Boussinesq)</td></tr>
    <tr><td><code>PARAM_T_END</code></td><td>real</td><td>&gt;0</td><td>Total simulation time</td></tr>
    <tr><td><code>PARAM_DT</code></td><td>real</td><td>&gt;0</td><td>Time step</td></tr>
    <tr><td><code>PARAM_CFL</code></td><td>real</td><td>0–1</td><td>CFL number (advection stability)</td></tr>
    <tr><td><code>PARAM_TIME_SCHEME</code></td><td>string</td><td>"RK4", "Euler", "RK2"</td><td>Time integration scheme</td></tr>
    <tr><td><code>PARAM_OUTPUT_EVERY</code></td><td>int</td><td>&ge;1</td><td>Save output every N steps</td></tr>
    <tr><td><code>PARAM_OUTPUT_FORMAT</code></td><td>string</td><td>"hdf5", "txt"</td><td>Output file format</td></tr>
    <tr><td><code>PARAM_OUTPUT_DIR</code></td><td>string</td><td>any valid folder name</td><td>Folder where results are saved</td></tr>
    <tr><td><code>PARAM_SAVE_VELOCITY</code></td><td>bool</td><td>TRUE / FALSE</td><td>Save velocity field</td></tr>
    <tr><td><code>PARAM_SAVE_PRESSURE</code></td><td>bool</td><td>TRUE / FALSE</td><td>Save pressure field</td></tr>
    <tr><td><code>PARAM_BC_TOP</code></td><td>string</td><td>"Dirichlet", "Neumann"</td><td>Top boundary condition</td></tr>
    <tr><td><code>PARAM_BC_BOTTOM</code></td><td>string</td><td>"Dirichlet", "Neumann"</td><td>Bottom boundary condition</td></tr>
    <tr><td><code>PARAM_BC_LEFT</code></td><td>string</td><td>"Dirichlet", "Neumann"</td><td>Left boundary condition</td></tr>
    <tr><td><code>PARAM_BC_RIGHT</code></td><td>string</td><td>"Dirichlet", "Neumann"</td><td>Right boundary condition</td></tr>
    <tr><td><code>PARAM_SOLVER_TYPE</code></td><td>string</td><td>"explicit", "implicit"</td><td>Time-stepping solver type</td></tr>
    <tr><td><code>PARAM_MAX_ITERS</code></td><td>int</td><td>&gt;0</td><td>Maximum iterations for iterative solvers</td></tr>
    <tr><td><code>PARAM_TOL</code></td><td>real</td><td>&gt;0</td><td>Convergence tolerance for solvers</td></tr>
    <tr><td><code>PARAM_LINEAR_SOLVER</code></td><td>string</td><td>"CG", "GMRES", "BiCGStab"</td><td>Linear solver algorithm</td></tr>
    <tr><td><code>PARAM_PRECONDITIONER</code></td><td>string</td><td>"ILU", "Jacobi", "None"</td><td>Preconditioner for linear solver</td></tr>
    <tr><td><code>PARAM_NUM_SCHEME</code></td><td>string</td><td>"central", "upwind"</td><td>Spatial discretization scheme</td></tr>
    <tr><td><code>PARAM_ORDER</code></td><td>int</td><td>1,2</td><td>Order of spatial scheme</td></tr>
    <tr><td><code>PARAM_STABILIZATION</code></td><td>real</td><td>&ge;0</td><td>Numerical stabilization (artificial diffusion)</td></tr>
    <tr><td><code>PARAM_VERBOSE</code></td><td>bool</td><td>TRUE / FALSE</td><td>Display console messages</td></tr>
    <tr><td><code>PARAM_CHECK_BOUNDS</code></td><td>bool</td><td>TRUE / FALSE</td><td>Check indices and grid boundaries</td></tr>
  </tbody>
</table>
