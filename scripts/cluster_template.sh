#!/bin/bash
#SBATCH --job-name={{ job_name }}
#SBATCH --partition=general
#SBATCH --output={{ log_file | /dev/null }} # you can add .%a for array index
#SBATCH --error={{ log_file | /dev/null }}
#SBATCH --mem-per-cpu={{ memory | 4096 }}
#SBATCH --array=1-{{ n_jobs }}
#SBATCH --cpus-per-task={{ cores | 1 }}

source activate perturbverse
# or: source activate {{ conda | default_conda_env_name }}
# or: your environment activation command
# ulimit -v $(( 1024 * {{ memory | 4096 }} ))
CMQ_AUTH={{ auth }} R --no-save --no-restore -e 'clustermq:::worker("{{ master }}")'
