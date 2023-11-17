#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=cbc-kinetic-transfer-kappa-1p5
#SBATCH --time=12:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

# `e281-ypimcf` is the Plasma Physics Consortium project code
#SBATCH --account=e281-ypimcf
#SBATCH --partition=standard
#SBATCH --qos=standard

# Set GS2 directory and run directory
GS2_exec="/work/e281/e281/bc1264/Applications/tobias-code/bin/gs2"
run_dir="/work/e281/e281/bc1264/tobias-code-initial-runs/cbc-electrostatic-kappa1p5/"
input_file="kappa1p5.in"
screen_output="screen_output.txt"

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

# Propagate the cpus-per-task setting from script to srun commands
#    By default, Slurm does not propagate this setting from the sbatch
#    options to srun commands in the job script. If this is not done,
#    process/thread pinning may be incorrect leading to poor performance
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Launch the parallel job
#   Using 256 MPI processes and 128 MPI processes per node
#   srun picks up the distribution from the sbatch options
srun --distribution=block:block --hint=nomultithread $GS2_exec $run_dir$input_file | tee $run_dir$screen_output