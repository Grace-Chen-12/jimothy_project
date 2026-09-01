#!/usr/bin/bash -l

#SBATCH --job-name=zarr_to_qupath
#SBATCH --partition=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/workdir/%u/logs/zarr_to_qupath_%j.out
#SBATCH --error=/workdir/%u/logs/zarr_to_qupath_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
#SBATCH --mail-user=your.email@example.com

mkdir -p /workdir/$USER/logs/

source /workdir/$USER/miniconda3/etc/profile.d/conda.sh # adjust to your conda install path
conda activate cellpose

NUCLEI_MASK="$1"
MICRONUCLEI_MASK="$2"
OUTPUT="$3"

ARGS=(--nuclei_mask "$NUCLEI_MASK")
[ -n "$MICRONUCLEI_MASK" ] && ARGS+=(--micronuclei_mask "$MICRONUCLEI_MASK")
[ -n "$OUTPUT" ] && ARGS+=(--output "$OUTPUT")

python -u zarr_to_qupath_annotations.py "${ARGS[@]}"
