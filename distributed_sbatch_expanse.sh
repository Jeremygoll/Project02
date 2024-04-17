#!/bin/bash
#SBATCH -p shared
#SBATCH -A mor101
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -D /home/colemans/project02

/home/colemans/project02/jeff_simulation-distributed -t $threads /home/colemans/project02/examples/$filename /dev/null