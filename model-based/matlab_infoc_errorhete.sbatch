#!/bin/bash 
#SBATCH --job-name=infocap_m
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=2
#SBATCH --output=jobid-%A_%a.out
#SBATCH --error=jobid-%A_%a.err
#SBATCH --array=1-100

#Load Matlab 2017a module

ml MATLAB

matlab -nojvm -nodisplay<<-EOF
for cond=1:2
infocapacity_hopf_errorhete(${SLURM_ARRAY_TASK_ID},cond);
end
EOF




#matlab -singleCompThread -nojvm -nodisplay < -r "run_hbif(${slurmArrayID},2,'grans');exit;" 
