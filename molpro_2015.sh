#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=1:00:00


# Controls the minimum/maximum number of nodes allocated to the job
#SBATCH -N 1

# Default resources are 1 core with 2.8GB of memory.

# Use more memory (4GB):
#SBATCH --mem=4G

# Specify a job name:
#SBATCH -J molpro

# Specify an output file
#SBATCH -o molpro-%j.out
#SBATCH -e molpro-%j.error

# Run a command
module load Molpro/2015_serial
molpro chd_reference.inp
