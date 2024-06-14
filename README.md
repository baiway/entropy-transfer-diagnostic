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

An example submission script is shown in [job.slurm](job.slurm). Ensure you change `run_dir` and `GS2_exec` before attempting to run, then submit the job with `sbatch job.slurm`.

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
Then create a virtual environment and install dependencies
```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
The [analysis](analysis/) folder contains the following scripts:
- [plot_poloidal_structure_Ts.py](analysis/plot_poloidal_structure_Ts.py) is a simple script that plots the poloidal structure of entropy transfer $T_s$.
- [movie-script.py](analysis/movie-script.py) produces a movie (.mp4) showing the evolution of several plots as the time integration period moves. The plots are: the $(k_x, k_y)$ spectrum of $\phi^2$, the time trace of $\phi^2$ and $T_s$, and the poloidal structure of $\phi^2$ and $T_s$. Detailed usage instructions can be found by entering `python3 movie-script.py --help`. **Beware:** this is not very carefully written and you'll run out of memory if the time integration window is too large. I recommend setting the window width (`T` in the code) to $\leq$ 200 (in GS2 units) if using inputs similar to those in [inputs.in](inputs.in).
- [dask-script.py](analysis/dask-script.py) I wrote this in an attempt to solve some of the memory issues I was encountering using the script above. It's currently very rough though and not finished. Would not recommend using yet.
- [dfn_plots.py](analysis/dfn_plots.py) uses the gs2ools community library to plot the distribution function as a function of $v_\parallel$ and $v_\perp$ for a particular $(k_x, k_y, \theta)$. Currently, the users specifies $(k_x, k_y)$ and the script produces velocity space plots for each $\theta$, then stiches them together into a movie.

## To-do:
- Add some notes explaining Eq. 16 and 17 in Ref. [2]
- In `src/diagnostics/gs2_diagnostics_new.f90`, I'm currently calling `init_diagnostics_entropy_transfer_3D` after `define_dims` (rather than with the other diagnostic initialisations). It would be better to move the dimension initilisation inside `define_dims` with a preprocessor flag (i.e. if using transfer diagnostic, define extra dimensions).

## References
[1] https://bitbucket.org/tobiasschuett/gs2-implement/

[2] M. Nakata, T.-H. Watanabe and H. Sugama, Nonlinear entropy transfer via zonal flows in gyrokinetic plasma turbulence, Physics of Plasmas 2012, https://aip.scitation.org/doi/full/10.1063/1.3675855
