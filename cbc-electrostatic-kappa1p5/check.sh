#!/bin/bash

# Check if a file argument is provided
if [ "$#" -ne 1 ]; then
    echo "Please provide a file argument. Usage: $0 <path-to-job.slurm>"
    exit 1
fi

# Extract variables from job.slurm using grep
GS2_exec=$(grep -oP '^GS2_exec="\K[^"]+' "$1")
run_dir=$(grep -oP '^run_dir="\K[^"]+' "$1")
input_file=$(grep -oP '^input_file="\K[^"]+' "$1")

# Check if the variables are set
unset_vars=()

[ -z "$GS2_exec" ] && unset_vars+=("GS2_exec")
[ -z "$run_dir" ] && unset_vars+=("run_dir")
[ -z "$input_file" ] && unset_vars+=("input_file")

if [ ${#unset_vars[@]} -ne 0 ]; then
    echo "The following variables are not set in $1: ${unset_vars[*]}"
    exit 1
fi

# 1) Check if run_dir exists
if [ -d "$run_dir" ]; then
    echo -e "\033[32mrun_dir exists: $run_dir\033[0m"
else
    echo -e "\033[31mrun_dir does not exist.\033[0m"
    exit 1
fi

# 2) Check if restart_dir exists within run_dir
restart_dir="${run_dir}restart_dir/"
if [ -d "$restart_dir" ]; then
    echo -e "\033[32mrestart_dir exists within run_dir.\033[0m"
else
    echo -e "\033[31mrestart_dir does not exist within run_dir.\033[0m"
fi

# 3) Check if input_file exists within run_dir
if [ -f "${run_dir}${input_file}" ]; then
    echo -e "\033[32minput_file exists within run_dir.\033[0m"
else
    echo -e "\033[31minput_file does not exist within run_dir.\033[0m"
fi

# 4) Check if GS2_exec is executable
if [ -x "$GS2_exec" ]; then
    echo -e "\033[32mGS2_exec is an executable: $GS2_exec\033[0m"
else
    echo -e "\033[31mGS2_exec is not an executable or does not exist.\033[0m"
fi
