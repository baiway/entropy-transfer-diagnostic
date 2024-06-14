# Entropy transfer workflow
Some brief documentation explaining how to use Tobias Schuett's entropy transfer diagnostic [1], based on Ref. [2].

## Getting started
### Compiling GS2
Clone GS2, switch to the `experimental/entropy-transfer-diagnostic` branch then update external submodules
```
git clone --recurse-submodules https://bitbucket.org/gyrokinetics/gs2.git
git checkout entropy-transfer-diagnostic
git submodule update --init --recursive
git submodule update --remote --recursive
git -C externals/utils checkout next
```

If running on Viking, load the following modules then compile:
```
module purge
module load gompi/2022b OpenMPI/4.1.4-GCC-12.2.0 # Requirements
module load netCDF-Fortran/4.6.0-gompi-2022b # For netcdf
module load FFTW/3.3.10-GCC-12.2.0 # For FFTW
module load OpenBLAS/0.3.21-GCC-12.2.0 # For Lapack/BLAS
module load Python/3.10.8-GCCcore-12.2.0 # For testing
module load compiler/GCC/7.3.0-2.30
module load mpi/OpenMPI/3.1.1-GCC-7.3.0-2.30
module load numlib/FFTW/3.3.8-gompi-2018b
module load data/netCDF-Fortran/4.4.4-foss-2018b
module load tools/git/2.19.1-GCCcore-7.3.0
make -I Makefiles GK_SYSTEM=viking -j 6
```

### Running the diagnostic
An example input file is shown in [input.in](input.in). I've just looked at the unshaped CBC so far with `tprim=4.0`. The settings for the diagnostic are:
```
&gs2_diagnostics_knobs
  ! ...
  write_entropy_transfer_3D = .true.
  entropy_transfer_symmetrised = .false.
/
```

The entropy transfer $T_s$ is written to the output NetCDF file. To access it using Xarray, for example, you can use:
```
import xarray as xr
ds = xr.open_dataset("output.nc.out", engine="netcdf4")
Ts = ds["entropy_transfer_3D"]
```

## Analysis
First clone this repository:
```
git clone https://github.com/baiway/entropy-transfer-diagnostic.git
```
The [analysis](analysis/) folder contains the following scripts:
- [plot_poloidal_structure_Ts.py](analysis/plot_poloidal_structure_Ts.py) is a simple script that plots the poloidal structure of entropy transfer $T_s$.
- Add movie script (+ tidy it up)
- Add Dask script (does it work?)


## To-do:
- In `src/diagnostics/gs2_diagnostics_new.f90`, I'm currently calling `init_diagnostics_entropy_transfer_3D` after `define_dims` (rather than with the other diagnostic initialisations). It would be better to move the dimension initilisation inside `define_dims` with a preprocessor flag (i.e. if using transfer diagnostic, define extra dimensions).

## References
[1] https://bitbucket.org/tobiasschuett/gs2-implement/
[2] M. Nakata, T.-H. Watanabe and H. Sugama, Nonlinear entropy transfer via zonal flows in gyrokinetic plasma turbulence, Physics of Plasmas 2012, https://aip.scitation.org/doi/full/10.1063/1.3675855
