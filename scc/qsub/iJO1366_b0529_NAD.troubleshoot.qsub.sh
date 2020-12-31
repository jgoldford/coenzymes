#!/bin/bash -l
# Set SCC project
#$ -P bioinfor

# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time.
#   The default time, also selected here, is 12 hours.  You can increase this up to 720:00:00 for single processor jobs but your job will take longer to start.
#$ -l h_rt=03:00:00

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

# Give job a name
#$ -N iJO1366_b0529_NAD_troubleshoot_job

# Combine output and error files into a single file
#$ -j y

module load miniconda
conda activate coenzymes

python /projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/coenzymes/lib/run_eFBA.py -m /projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/coenzymes/assets/models/stoich_thermo/iJO1366_b0529_NAD.xml -e /projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/coenzymes/assets/media.C180N93E7.txt -o /projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/coenzymes/scc/output/iJO1366_b0529_NAD.ts.txt