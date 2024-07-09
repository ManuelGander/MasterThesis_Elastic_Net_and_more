#!/bin/bash
#
#SBATCH --job-name=Moscot_job
#SBATCH --output=/home/icb/manuel.gander/do_k.txt
#SBATCH --error=/home/icb/manuel.gander/error_k.txt
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=10G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000


source /home/icb/manuel.gander/ott_env/bin/activate


i="$1"
j="$2"
k="$3"

python ps6.py $i $j $k





