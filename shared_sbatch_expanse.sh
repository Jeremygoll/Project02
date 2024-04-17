#!/bin/bash
#SBATCH -A mor101
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -D /home/colemans/project02

/home/colemans/project02/simulation-shared -t $threads /home/colemans/project02/examples/$filename /dev/null