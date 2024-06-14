#!/usr/bin/env bash
# SBATCH --job-name=entropy-transfer    # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=00-48:00:00              # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=240                    # Number of MPI tasks to request
#SBATCH --cpus-per-task=1               # Number of CPU cores per MPI task
#SBATCH --mem-per-cpu=5200              # Default mem-per-cpu (default units are megabytes)
#SBATCH --account=pet-gspt-2019         # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

# Abort if any command fails
set -e

# Purge any previously loaded modules
module purge

# Load modules
module load gompi/2022b OpenMPI/4.1.4-GCC-12.2.0 # Requirements
module load netCDF-Fortran/4.6.0-gompi-2022b     # For netcdf
module load FFTW/3.3.10-GCC-12.2.0               # For FFTW
module load OpenBLAS/0.3.21-GCC-12.2.0           # For Lapack/BLAS
module load Python/3.10.8-GCCcore-12.2.0         # For testing

# Extract run_dir from the path to job.slurm
run_dir="/users/bc1264/scratch/cbc-distfn-transfer"

# Set GS2 directory and run directory
GS2_exec="/users/bc1264/scratch/gs2-transfer-neasy-f/bin/gs2"
input_file="input.in"
screen_output="screen_output.txt"

# Commands to run
srun $GS2_exec $run_dir/$input_file | tee $run_dir/$screen_output
