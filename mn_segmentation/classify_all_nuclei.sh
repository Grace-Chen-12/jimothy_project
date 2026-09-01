#!/usr/bin/bash -l

#SBATCH --job-name=classify_all_nuclei
#SBATCH --partition=regular
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/workdir/%u/logs/classify_nuclei_%j.out
#SBATCH --error=/workdir/%u/logs/classify_nuclei_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
#SBATCH --mail-user=your.email@example.com

# change the sbatch directives to match your cluster's configuration

mkdir -p /workdir/$USER/logs/

# change conda environment name to match your setup
source /workdir/$USER/miniconda3/etc/profile.d/conda.sh 
conda activate cellpose

IMAGE="$1"
OUTPUT="$2"
CHANNEL_MIN="$3"
CHANNEL_MAX="$4"
PRETRAINED_MODEL="${5:-cpsam_v2}"

ARGS=(--image "$IMAGE")
[ -n "$OUTPUT" ] && ARGS+=(--output "$OUTPUT")
[ -n "$CHANNEL_MIN" ] && ARGS+=(--channel_min "$CHANNEL_MIN")
[ -n "$CHANNEL_MAX" ] && ARGS+=(--channel_max "$CHANNEL_MAX")
[ -n "$PRETRAINED_MODEL" ] && ARGS+=(--pretrained_model "$PRETRAINED_MODEL")

python -u classify_all_nuclei.py "${ARGS[@]}"