#!/bin/bash
#
#SBATCH --job-name=Moscot_job
#SBATCH --output=/home/icb/manuel.gander/do_k.txt
#SBATCH --error=/home/icb/manuel.gander/error_k.txt
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --nice=10000


source /home/icb/manuel.gander/ott_env/bin/activate


i="$1"
j="$2"
k="$3"
l="$4"

python ps0.py $i $j $k $l