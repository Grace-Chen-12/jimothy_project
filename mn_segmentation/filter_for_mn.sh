#!/usr/bin/bash -l

#SBATCH --job-name=filter_for_mn
#SBATCH --partition=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/workdir/%u/logs/filter_for_mn_%j.out
#SBATCH --error=/workdir/%u/logs/filter_for_mn_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
#SBATCH --mail-user=your.email@example.com

mkdir -p /workdir/$USER/logs/

# adjust sbatch parameters and conda environment as needed
source /workdir/$USER/miniconda3/etc/profile.d/conda.sh 
conda activate cellpose

NUCLEI_MASK="$1"
CELLS_ZARR="$2"
OUTPUT="${3:-}"
POST_IF_ZARR="${4:-}"

ARGS=(--nuclei_mask "$NUCLEI_MASK" --cells_zarr "$CELLS_ZARR")
[ -n "$OUTPUT" ] && ARGS+=(--output "$OUTPUT")
[ -n "$POST_IF_ZARR" ] && ARGS+=(--post_if_zarr "$POST_IF_ZARR")

python -u filter_for_mn.py "${ARGS[@]}"
