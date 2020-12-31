import glob
import subprocess as sb

def make_qsub(job_name,eFBAcall,modelFile,mediaFile,outFile):
    out_string = f"""#!/bin/bash -l
# Set SCC project
#$ -P bioinfor

# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time.
#   The default time, also selected here, is 12 hours.  You can increase this up to 720:00:00 for single processor jobs but your job will take longer to start.
#$ -l h_rt=12:00:00

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

# Give job a name
#$ -N {job_name}

# Combine output and error files into a single file
#$ -j y

module load miniconda
conda activate coenzymes
module load gurobi

"""
    python_call = 'python ' + eFBAcall + ' -m ' + modelFile + ' -e ' + mediaFile + ' -o ' + outFile;
    out_string = out_string + python_call
    return out_string

#root_dir = '/Users/Joshua.Goldford/Documents/GitHub/coenzymes/'
root_dir = '/projectnb/bioinfor/SEGRE/goldford/CoenzymeSpecificity/coenzymes/'
python_efba_path = root_dir + 'lib/run_eFBA.py'
media_file = root_dir + 'assets/media.C180N93E7.txt'


# make and submit jobs for stoich only models

models = glob.glob(root_dir + 'assets/models_dg_range_2/stoich/*.xml')

for model in models:
    name = model.split('/')[-1].split('.')[0]
    name = name + '_S'
    outFile = root_dir + 'scc/output_bounds2_gurobi/' + name + '.txt'
    qsub = make_qsub(name,python_efba_path,model,media_file,outFile)
    sh_file = root_dir + 'scc/qsub_bounds2_gurobi/' + name + '.qsub'
    
    with open(sh_file,'w') as f:
        f.write(qsub)

    # call qsub file to submit job
    sb.call('qsub ' + sh_file,shell=True)
    
# make and submit jobs for stoich + thermo models

models = glob.glob(root_dir + 'assets/models_dg_range_2/stoich_thermo/*.xml')

for model in models:
    name = model.split('/')[-1].split('.')[0]
    name = name + '_ST'
    outFile = root_dir + 'scc/output_bounds2_gurobi/' + name + '.txt'
    qsub = make_qsub(name,python_efba_path,model,media_file,outFile)
    sh_file = root_dir + 'scc/qsub_bounds2_gurobi/' + name + '.qsub'
    
    with open(sh_file,'w') as f:
        f.write(qsub)

    # call qsub file to submit job
    sb.call('qsub ' + sh_file,shell=True)