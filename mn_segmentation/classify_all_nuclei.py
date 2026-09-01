'''
Goal: Segment all nuclei (PN and MN) using Cellpose.


Input file: ch0000_dapi.ome.tif file from morphology_focus folder of Xenium outputs
Output file: Zarr file with Cellpose segmentation
'''
import numpy as np
import tifffile
import zarr
from cellpose import models
from skimage.measure import regionprops
import argparse
import os
from pathlib import Path

def open_image(path):
    '''Open zarr view of tiff store'''
    store = tifffile.imread(path, aszarr=True, is_ome=False) # is_ome=False stops tiffile from re-merging per-channel series into a single array
    root = zarr.open(store, mode='r')
    if isinstance(root, zarr.hierarchy.Group):
        return root['0'] # extracts 0th (DAPI) channel
    return root

def auto_channel_range(img, low_percentile=1.0, high_percentile=99, sample_step=8):
    '''Estimates a DAPI display range (floor, ceiling) from a strided subsample image'''
    sample = np.asarray(img[::sample_step, ::sample_step])
    lo, hi = np.percentile(sample, (low_percentile, high_percentile))
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)

def make_tiles(height, width, tile_size, overlap):
    """Yield (y0, y1, x0, x1) bounds covering (height, width) in tile_size
    chunks with the given overlap between neighboring tiles."""
    step = tile_size - overlap
    for y0 in range(0, height, step):
        y1 = min(y0 + tile_size, height)
        for x0 in range(0, width, step):
            x1 = min(x0 + tile_size, width)
            yield y0, y1, x0, x1
            if x1 == width:
                break
        if y1 == height:
            break

def build_cache(img, cache_path, tile_size, source_path):
    """Decode image once and store in a zarr cache"""
    cache_path = Path(cache_path)
    source_mtime = os.path.getmtime(source_path)

    # If the cache already exists and from same image source file, just return it
    if cache_path.exists():
        cache = zarr.open(str(cache_path), mode='r')
        if cache.attrs.get('source_mtime') == source_mtime:
            print(f"Using existing cache at {cache_path}")
            return cache
        print(f"Cache at {cache_path} is stale (source file changed). Rebuilding...")
    
    # Otherwise, build the cache from scratch
    height, width = img.shape[0], img.shape[1]
    cache = zarr.open(
        str(cache_path), mode='w', shape=(height, width), chunks=(tile_size, tile_size), dtype=img.dtype,
    )
    cache.attrs['source_mtime'] = source_mtime

    print(f"Building cache of input image at {cache_path}...")

    # Decompresses TIFF data once so we don't have to redo decode for every tile
    for y0 in range(0, height, tile_size):
        y1 = min(y0 + tile_size, height)
        cache[y0:y1, :] = np.asarray(img[y0:y1, :])
    print("Done building cache.")

    # Return read-only view of the cache so we don't accidently overwrite
    return zarr.open(str(cache_path), mode='r')

def remove_spurious_annotations(label_image, intensity_image, min_mean_intensity):
    """Zero out labels whose mean DAPI intensity falls below min_mean_intensity."""
    filtered = label_image.copy()
    n_removed = 0
    for prop in regionprops(label_image, intensity_image=intensity_image):
        if prop.mean_intensity < min_mean_intensity:
            rows, cols = prop.coords[:, 0], prop.coords[:, 1]
            filtered[rows, cols] = 0
            n_removed += 1
    return filtered, n_removed

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--image', type=str, required=True, help='Path to the input DAPI image (ch0000_dapi.ome.tif)')
    parser.add_argument('--output', type=str, default=None, help='Path to write the nuclei label mask (ends with .zarr)')
    parser.add_argument('--channel_min', type=float, default=None,
                        help='DAPI normalization floor, about 1.0 is generally good')
    parser.add_argument('--channel_max', type=float, default=None,
                        help='DAPI normalization ceiling, about 400.0 is generally good')
    parser.add_argument('--pretrained_model', type=str, default='cpsam_v2',
                        help='Cellpose model name or path to a custom model')
    return parser.parse_args()


def main():
    print("Running")

    args = parse_args()

    IMAGE = args.image
    OUTPUT = args.output or str(Path(IMAGE).with_suffix('')) + '_nuclei_mask.zarr'

    img = open_image(IMAGE)
    assert img.ndim == 2, f"Expected 2D image after processing, got {img.ndim}D image instead"
    HEIGHT, WIDTH = img.shape # Y, X
    print(f"Image shape (Y, X): {img.shape}")

    # Use user-specified channel min/max if provided, otherwise auto-compute
    if args.channel_min is not None and args.channel_max is not None:
        CHANNEL_MIN = args.channel_min
        CHANNEL_MAX = args.channel_max
    else:
        CHANNEL_MIN, CHANNEL_MAX = auto_channel_range(img)
        print(f"channel_min/channel_max not specified; using auto-computed range: "
        f"{CHANNEL_MIN:.1f} / {CHANNEL_MAX:.1f}")

    # Hardcoded parameters
    TILE_SIZE = 4096
    OVERLAP = 256
    MIN_MEAN_DAPI_INTENSITY = 25.0 # labels below this threshold are treated as spurious

    # Cellpose parameters
    FLOW_THRESHOLD = 0.5 # default = 0.4
    CELLPROB_THRESHOLD = 0.5 # default = 0.0
    MIN_SIZE = 4 # default = 15
    PATH_TO_MODEL = args.pretrained_model

    tiles = list(make_tiles(HEIGHT, WIDTH, TILE_SIZE, OVERLAP))
    print(f"Generated {len(tiles)} tiles of size {TILE_SIZE}px with {OVERLAP}px overlap")

    CACHE_PATH = str(Path(IMAGE).with_suffix('')) + '_cache.zarr'
    dapi = build_cache(img, CACHE_PATH, TILE_SIZE, IMAGE)


    model = models.CellposeModel(gpu=True, pretrained_model=PATH_TO_MODEL)

    mask_store = zarr.open(
        OUTPUT, mode='w', shape=(HEIGHT, WIDTH), chunks=(TILE_SIZE, TILE_SIZE), dtype=np.uint32
    )


    label_id = 0 # global label ID counter across all tiles
    for i, (y0, y1, x0, x1) in enumerate(tiles):
        # Per tile normalization to [0, 1]
        LOW, HIGH = CHANNEL_MIN, CHANNEL_MAX
        dapi_tile = np.asarray(dapi[y0:y1, x0:x1])
        tile = dapi_tile.astype(np.float32)
        tile = np.clip((tile - LOW) / (HIGH - LOW), 0.0, 1.0)

        masks, flows, styles = model.eval(
            tile,
            diameter=None,
            normalize=False, # already normalized earlier
            flow_threshold=FLOW_THRESHOLD,
            cellprob_threshold=CELLPROB_THRESHOLD,
            min_size=MIN_SIZE,
        )

        # This tile "owns" the half-overlap core region below.
        core_y0 = OVERLAP // 2 if y0 > 0 else 0
        core_x0 = OVERLAP // 2 if x0 > 0 else 0
        core_y1 = (y1 - y0) - (OVERLAP // 2 if y1 < HEIGHT else 0)
        core_x1 = (x1 - x0) - (OVERLAP // 2 if x1 < WIDTH else 0)

        # Build the relabeled result for this tile in a plain numpy array 
        n_kept = 0 # per-tile count for reporting
        n_spurious = 0

        tile_result = np.zeros(masks.shape, dtype=np.uint32)
        for prop in regionprops(masks):
            # to account for nuclei spanning multiple tiles, the tile that
            # "owns" the nucleus is based on the nucleus centroid
            cy, cx = prop.centroid
            if not (core_y0 <= cy < core_y1 and core_x0 <= cx < core_x1):
                continue  # centroid belongs to a neighboring tile
            n_kept += 1
            label_id += 1
            rows, cols = prop.coords[:, 0], prop.coords[:, 1]
            tile_result[rows, cols] = label_id

        if n_kept > 0:
            tile_result, n_spurious = remove_spurious_annotations(tile_result, dapi_tile, MIN_MEAN_DAPI_INTENSITY)

            # merge rather than overwrite, don't lose existing labels from neighboring tiles
            existing = np.asarray(mask_store[y0:y1, x0:x1])
            mask_store[y0:y1, x0:x1] = np.where(existing > 0, existing, tile_result)

        n_found = int(masks.max())

        print(f"Tile {i + 1}/{len(tiles)} ({y0}:{y1}, {x0}:{x1}): {n_kept} nuclei kept (of {n_found} found), "
            f"{n_spurious} spurious annotations removed (mean DAPI < {MIN_MEAN_DAPI_INTENSITY})")


    print(f"Nuclei mask saved to {OUTPUT}")


if __name__ == "__main__":
    main()