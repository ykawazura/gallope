# <span style="font-variant: small-caps;">Gallope</span>

![Language](https://img.shields.io/badge/language-Fortran-purple)
![Acceleration](https://img.shields.io/badge/Acceleration-OpenACC%20%7C%20cuFFTMp-green)
![License](https://img.shields.io/badge/license-MIT-blue)

**<span style="font-variant: small-caps;">Gallope</span>** is a multi-node, multi-GPU pseudo-spectral code designed for high-performance simulations of Magnetohydrodynamics (MHD). It is developed from the CPU code **[<span style="font-variant: small-caps;">Calliope</span>](https://github.com/ykawazura/calliope)** with minimal modifications, accelerated using **OpenACC** and **NVIDIA cuFFTMp**.

---

## Highlights

- Multi-node multi-GPU FFT via NVIDIA cuFFTMp
- OpenACC implementation (keeps CPU/GPU code divergence small)
- Designed for large 3D grids and strong scaling on GPU clusters
- Minimal changes from <span style="font-variant: small-caps;">Calliope</span> for maintainability

---

## Quick start
```bash
git clone https://github.com/ykawazura/gallope.git
cd gallope

# 1. Configure Makefile.in and arch/<arch>.in (see Build section)
# 2. Compile
make

# 3. Prepare input
cp "input_example/MHD_INCOMP(forcing).in" ./gallope.in

# 4. Run
mpirun -np <NRANK> ./gallope
```

---

## Models

Set `MODEL` in `Makefile.in`:

| `MODEL`        | Description        |
|----------------|--------------------|
| `MHD_INCOMP`   | Incompressible MHD |
| `RMHD`         | Reduced MHD        |


---

## Requirements

### Compilers & Toolchain
- NVIDIA HPC SDK (NVHPC)
    - Must support OpenACC.
- NVIDIA cuFFTMp (Required for multi-node FFT)
- MPI Implementation

### I/O Libraries
- netCDF (C library)
- netCDF-Fortran

> **Note:** <span style="font-variant: small-caps;">Gallope</span> links both `-lnetcdf` and `-lnetcdff`. Depending on your specific build, you may also need to link HDF5.

---

## Build Instructions

### 1. Configure `Makefile.in`
Set the target environment (`arch`) and the physics model (`MODEL`).

```make
arch = miyabi           # Loads configuration from arch/miyabi.in
MODEL = MHD_INCOMP      # Options: MHD_INCOMP, RMHD
````

### 2. Configure `arch/<arch>.in`

Create or edit a file in the `arch/` directory (e.g., `arch/miyabi.in`). You must define the compiler, library paths, and flags.

Example (`arch/miyabi.in`):

```make
F90 = mpif90

# NetCDF Paths
NETCDF_INC = -I$(NETCDF_DIR)/include -I$(NETCDF_FORTRAN_DIR)/include
NETCDF_LIB = -L$(NETCDF_DIR)/lib -L$(NETCDF_FORTRAN_DIR)/lib -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lz

# Compilation Flags (NVHPC + OpenACC + cuFFTMp)
FLAGS  = -O3 -Mfree -fast -Mextend -Mpreprocess -Minform=warn -Minfo=accel -cuda -cudalib=cufftmp
LINKER = -L${NVHPC_ROOT}/compilers/lib -lnvhpcwrapcufft -L${NVHPC_ROOT}/math_libs/lib64/ -lcufftMp -lmpi
```


Adjust `NETCDF_DIR` and `NETCDF_FORTRAN_DIR` to match your environment. Additional libraries (e.g., `-ldl`, `-lm`) may be needed.

### 3. Compile

```bash
make
```

---

## Execution

### 1. Prepare Input File

<span style="font-variant: small-caps;">Gallope</span> reads configuration from a file named **`gallope.in`** in the current directory. Copy an example to get started:


```bash
# Example: Incompressible MHD with forcing
cp input_example/MHD_INCOMP\(forcing\).in ./gallope.in
```

### 2. Run

Execute with MPI. 

```bash
mpirun -np <NRANK> ./gallope
```

### 3. Example Standard Output

```Plaintext
 --------------------------------------------------------------------
|                                                                    |
|     Gallope started!                                               |
|                                     executed on 2026- 1- 7 13:33   |
 --------------------------------------------------------------------

 ==========================================
 MPI processes:      16
 Nodes:              16
 GPUs per node:       1
 ==========================================
 This run will be terminated when the wall time exceeds  9.50 hours.

 nlx  =   1024,    nly  =   1024,    nlz  =    512

Solving MHD_INCOMP
Restart
shear is on; tremap = 1.333
----------------------------------------------------------------------------
Initialization done
Starting the main loop...
step =        551/    100000,  time =      36.79988937
...
step =      26405/    100000,  time =      20.29980469

!-----------------------------------------------------------------------!
!  Wall time has exceeded the set value. Terminating the simulation...  !
!-----------------------------------------------------------------------!
Finished the main loop :)

Initialization                 0.172 min     0.1 %
Advance steps                263.590 min    97.5 %
   nonlinear terms           236.067 min    87.3 %
   FFT                       220.391 min    81.5 %
   IO                          0.038 min     0.0 %
Save restart                   1.634 min     0.6 %

total from timer is:     270.398 min
 ---  Program completed! ^_^ ---
```

---

## Output Files

|**File/Directory**|**Description**|
|---|---|
|**`gallope.out.nc`**|netCDF file containing global diagnostics (energies) and 1D spectra.|
|**`out2d/`**|Time evolution of 2D slices (cross-sections) of fields.|
|**`out3d/`**|Time evolution of full 3D field data.|
|**`restart/`**|Snapshots generated at the end of the run (or periodically).|

To restart a simulation, set `init_type = 'restart'` in your `gallope.in` file.

---

## Diagnostics & Visualization

Python visualization scripts are provided in the `diagnostics/` directory, organized by model.

**Location:** `diagnostics/<MODEL>/` (e.g., `diagnostics/MHD_INCOMP/`)

| Script              | Description                        |
|---------------------|------------------------------------|
| **`plot_energy.py`**    | Time evolution of energies         |
| **`plot_kspectrum.py`** | 1D energy spectra                  |
| **`plot_fields.py`**    | 2D field slices    
    

---

---

## Reproducing the paper

The measurements behind *"Gallope: A Multi-GPU pseudospectral solver for
magnetohydrodynamic models using the cuFFTMp library"* (COMPHY-D-26-00110) are archived, not summarised.

| what | where |
|---|---|
| Section 4, Figures 1-5 — FFT-level and application-level scaling on Miyabi-G and TSUBAME 4.0 | **[`benchmarks/`](benchmarks/)** in this repository: raw stdout of every run, the inputs and job scripts that produced it, the campaign drivers, and the plotting scripts |
| Section 5, Figures 6-8 and Table 2 — deterministic <span style="font-variant: small-caps;">Gallope</span>/<span style="font-variant: small-caps;">Calliope</span> comparison and the nonlinear MRI run | **[ykawazura/gallope-calliope_comparison](https://github.com/ykawazura/gallope-calliope_comparison)**, where each run sits beside its <span style="font-variant: small-caps;">Calliope</span> counterpart |

`benchmarks/README.md` carries the figure-by-figure index, the correspondence
to the benchmarking protocol of Section 4.1, and how to re-plot each figure.

The code those measurements were taken with is this branch at tag
[`v1.0`](https://github.com/ykawazura/gallope/releases/tag/v1.0).  Cite the tag
rather than the tip of `main`: later commits will move on from the paper.

## Repository Structure

```Plaintext
.
├── src/                # Source code
├── arch/               # Machine-dependent compile settings
│   └── miyabi.in
│   └── ...
├── input_example/      # Sample input files
│   └── MHD_INCOMP(forcing).in
│   └── ...
├── diagnostics/        # Analysis scripts 
│   ├── MHD_INCOMP/
│   └── RMHD/
│   └── ...
├── benchmarks/         # Raw timings, inputs and plotting scripts of the paper's Figs. 1-5
│   ├── gallope/                # application-level scaling
│   ├── p3dfft_vs_cufftmp/      # FFT-level scaling
│   └── README.md
├── tests/              # Reference runs with their recorded stdout
├── Makefile.in         # Main build 
└── README.md
```

---

<!--
## Citation

If you use <span style="font-variant: small-caps;">Gallope</span> in academic work, please cite:

1. **[Paper Title / Journal]** (TODO: Add Link/DOI)
    
2. **Code Archive:** (TODO: Add Zenodo DOI)
    

Please also consider citing the original **[<span style="font-variant: small-caps;">Calliope</span> paper](https://doi.org/10.3847/1538-4357/ac4f63)** paper if applicable.

---
-->

## Contact

**Maintainer:** Yohei Kawazura (kawazura@utsunomiya-u.ac.jp)

---

## License

[MIT License](https://www.google.com/search?q=LICENSE)
