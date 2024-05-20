#!/bin/bash
#SBATCH --job-name=seurat-anchors
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=100GB
#SBATCH --array=0-26

__conda_setup="$('/camp/home/coultoa/homehd/anaconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/camp/home/coultoa/homehd/anaconda/etc/profile.d/conda.sh" ]; then
        . "/camp/home/coultoa/homehd/anaconda/etc/profile.d/conda.sh"
    else
        export PATH="/camp/home/coultoa/homehd/anaconda/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate renv
Rscript ~/work/ucl/scripts/ucl.projects/macrophage/find.markers2.R ${SLURM_ARRAY_TASK_ID} 
