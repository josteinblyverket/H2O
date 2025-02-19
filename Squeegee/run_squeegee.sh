#!/usr/bin/bash
#SBATCH --job-name=squeegee
#SBATCH --qos=nf
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
##SBATCH --mem=6GB
#SBATCH --account=nores
##SBATCH --output=g.%j.out
##SBATCH --error=/ec/res4/scratch/sbjb/Projects/CARRA2/logfiles/copy_forcing.%j.out

module load python3/3.11.8-01
cd /ec/res4/hpcperm/sbjb/github/H2O/Squeegee/Squeegee

python3 cerise_runoff.py
