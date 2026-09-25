#!/bin/bash

# Copy this file to a new job script and update the marked values.
# Replace <repo-root> with the absolute path to your repository checkout.

#SBATCH --job-name="ExampleJob"
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --account=<project-account>
#SBATCH --chdir=<repo-root>
#SBATCH --output=<repo-root>/slurm_jobs/%x-%j.out
#SBATCH --error=<repo-root>/slurm_jobs/%x-%j.err

# Replace this with the matching run_slurm entrypoint for your job.
bash run_slurm/your_script.sh