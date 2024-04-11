#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 15                      # Number of tasks (really number of CPU Cores/task). Don't specify more than 16 unless approved by the system admin
#SBATCH --array=28-42    #change this to match up with the number of parallel computing jobs; also generate a config file 
#SBATCH --mail-type=END           # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE 
#SBATCH --mail-user=kexindon@mit.edu           # Email to which notifications will be sent. Equivalent to the -M option in SGE. You must replace [] with your email address.
#############################################
srun --pty bash

cd /net/bmc-lab2/data/lab/sanchezrivera/kexindong

module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit

# conda env create -n new_h2m_env -f ./h2m_env.yml

#----then in the future to actiavte-------
# module load miniconda3/v4
# source /home/software/conda/miniconda3/bin/condainit
conda activate new_h2m_env

cd /net/bmc-lab2/data/lab/sanchezrivera/kexindong/h2m-aacr-march-25

#access the config file
config=./h2m_config.txt

FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
folder_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

# python3 crispresso_analysis_aggregation_42nt.py ${FILE} ${folder_name}

#test 
# echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${FILE} and the output folder is ${folder_name}" >> output.txt
# pip3 install --force-reinstall h2m-0.1.22-py3-none-any.whl
python3 task_aacr.py ${FILE} ${folder_name}