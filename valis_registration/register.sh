#!/usr/bin/bash -l

#SBATCH --job-name=valis_register
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=350G
#SBATCH --time=24:00:00
#SBATCH --output=/workdir/%u/logs/register_%j.out
#SBATCH --error=/workdir/%u/logs/register_%j.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
#SBATCH --mail-user=your.email@example.com

# change --mail-user in the sbatch directives to your own email


# Directory the code lives in
CODE_DIR="${CODE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"

# strip the /local prefix from /local/workdir to get a path the container can see.
CODE_DIR="${CODE_DIR#/local}"

# Host directory made visible inside the container
MOUNT="/workdir/$USER"

# Mount both spellings so either form works for the paths you pass in.
MOUNTS=(-v "$MOUNT:$MOUNT")
[ -d "/local$MOUNT" ] && MOUNTS+=(-v "/local$MOUNT:/local$MOUNT")

# Run VALIS inside Docker.
# --no-deps is required here as imagecodecs depends on numpy and without this,
# pip would upgrade the image's numpy past the numpy<2.0.0 that valis-wsi pins.
docker1 run --rm --memory=350g \
  -e PYTHONUNBUFFERED=1 \
  "${MOUNTS[@]}" \
  cdgatenbee/valis-wsi:1.2.0 \
  bash -c 'pip install --quiet --no-cache-dir --no-deps imagecodecs && exec python3 "$0" "$@"' \
  "$CODE_DIR/register.py" "$@"
