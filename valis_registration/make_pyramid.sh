#!/bin/bash
#SBATCH --job-name=ome_pyramid
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=8:00:00
#SBATCH --output=/workdir/%u/logs/pyramid_%j.out
#SBATCH --error=/workdir/%u/logs/pyramid_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
#SBATCH --mail-user=your.email@example.com

# change --mail-user in the sbatch directives to your own email

set -euo pipefail
mkdir -p /workdir/$USER/logs

python make_pyramid.py "$1"