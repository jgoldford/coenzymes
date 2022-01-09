#!/bin/bash -l
# Set SCC project
#$ -P bioinfor

# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time.
#   The default time, also selected here, is 12 hours.  You can increase this up to 720:00:00 for single processor jobs but your job will take longer to start.
#$ -l h_rt=12:00:00

# Give job a name
#$ -N tfa_v1

# Combine output and error files into a single file
#$ -j y

# Add PE specifications
#$ -pe omp 4


module load miniconda
conda activate coenzymes
module load gurobi

python /projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/pytfa/tfa_sc_results/run_tfa_wt_sc_model.py
