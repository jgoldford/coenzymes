#!/bin/bash -l

# Set SCC project
#$ -P bioinfor

# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time.
#   The default time, also selected here, is 12 hours.  You can increase this up to 720:00:00 for single processor jobs but your job will take longer to start.
#$ -l h_rt=01:00:00

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

# Give job a name
#$ -N test_eFBA_10envs

# Combine output and error files into a single file
#$ -j y

module load miniconda
conda activate coenzymes

# run the python script
python /projectnb/bioinfor/goldford/tests/run_eFBA.py -m /projectnb/bioinfor/goldford/tests/iJO1366.xml -e /projectnb/bioinfor/goldford/tests/media.C180N93E7.sample1000.txt -o /projectnb/bioinfor/goldford/tests/test1000.txt 
