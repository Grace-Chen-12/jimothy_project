"""
Perform IF-on-IF or HE-on-IF WSI registration using VALIS.

Submit with register.sh
"""

import argparse
import os
import time
from pathlib import Path
import numpy as np

import io_utils
import stains
import transforms
import valis_runner

def tag_filename(name, *tags):
    """
    Insert "tags" before the file extension(s).
    Example:
        tag_filename("image.ome.tif", "staged", "rot90")
        -> "image_staged_rot90.ome.tif"
    """
    p = Path(name)
    # ext holds concatenated suffixes (e.g., ".ome.tif")
    ext = "".join(p.suffixes)
    # stem holds filename w/out any suffixes
    stem = p.name[: -len(ext)] if ext else p.name
    return f"{stem}_{'_'.join(tags)}{ext}"

def parse_args():
    parser = argparse.ArgumentParser(description="WSI registration")
    parser.add_argument("--moving_image", type=Path, required=True, help="Path to the moving image")
    parser.add_argument("--fixed_image", type=Path, required=True, help="Path to the fixed image (reference, always IF)")
    parser.add_argument("--output_dir", type=Path, required=True, help="Directory to save the registered images and any intermediate files")
    parser.add_argument("--moving_modality", choices=["HE", "IF"], required=True, help="Modality of the moving image. Reference is always IF.")
    parser.add_argument("--transform", type=str, default="identity", help= "one of: identity, rot90, rot180, rot270, flipud, fliplr, transpose, transverse (default: identity)")
    return parser.parse_args()

def main():
    args = parse_args()
    moving_image = args.moving_image
    fixed_image = args.fixed_image
    output_dir = args.output_dir
    moving_modality = args.moving_modality
    transform = args.transform

    compression = "LZW" # I've found LZW compression to work with VALIS

    print("\n--- Inputs ---")
    print(f"moving image   : {moving_image}")
    print(f"fixed image    : {fixed_image}")
    print(f"output dir     : {output_dir}")
    print(f"moving modality: {moving_modality} (reference is always IF)")
    print(f"transform      : {transform}")
    print(f"compression    : {compression}")

    # Read the moving image together with its OME-XML metadata
    print("\n--- Reading moving image ---")
    start = time.time()
    moving_arr, moving_xml, axes = io_utils.read_ome(moving_image)
    y_axis_idx, x_axis_idx = io_utils.spatial_axes(axes)
    print(f"shape {moving_arr.shape}, dtype {moving_arr.dtype}, axes '{axes}' (Y={y_axis_idx}, X={x_axis_idx})")
    if moving_xml is None:
        print("WARNING: no OME-XML metadata found, staged image will have none either")
    print(f"read in {(time.time() - start)/60:.2f} minutes")

    # Ensure that for H&E images, the channel axis is last
    if moving_modality == "HE":
        c_axis_idx = io_utils.channel_axis(axes)
        if c_axis_idx is None:
            raise ValueError(f"moving modality is HE but axes '{axes}' has no channel axis")
        if c_axis_idx != moving_arr.ndim - 1:
            print(f"\n--- Moving channel axis {c_axis_idx} last ---")
            moving_arr = np.moveaxis(moving_arr, c_axis_idx, -1)
            axes = axes.replace(axes[c_axis_idx], "") + axes[c_axis_idx]
            y_axis_idx, x_axis_idx = io_utils.spatial_axes(axes)
            # The samples are interleaved now whatever the source file said
            moving_xml = io_utils.set_interleaved(moving_xml, True)
            print(f"shape {moving_arr.shape}, axes '{axes}' (Y={y_axis_idx}, X={x_axis_idx})")

    # Transform pixels and metadata
    print(f"\n--- Applying transform '{transform}' ---")
    moving_arr = transforms.apply_transform(moving_arr, transform, y_axis_idx, x_axis_idx)
    print(f"shape is now {moving_arr.shape}")
    if transforms.swaps_pixels(transform):
        moving_xml = io_utils.swap_xy_metadata(moving_xml)

    # Make a directory for the staged images in the output directory
    staged_dir = output_dir / "staged"
    os.makedirs(staged_dir, exist_ok=True)
    # Build filenames for staged images
    moving_name = tag_filename(moving_image.name, transform, "staged")
    reference_name = tag_filename(fixed_image.name, "staged")

    # Stage the transformed RGB, 
    # then deconvolve the hematoxylin channel of H&E
    rgb_staged_path = None
    if moving_modality == "HE":
        rgb_dir = output_dir / "staged_rgb"
        os.makedirs(rgb_dir, exist_ok=True)
        rgb_staged_path = rgb_dir / tag_filename(moving_image.name, transform, "rgb", "staged")
        print(f"\n--- Staging transformed RGB ---")
        print(f"rgb: {rgb_staged_path}")
        start = time.time()
        io_utils.write_ome_tif(str(rgb_staged_path), moving_arr, moving_xml, compression=compression)
        print(f"staged in {(time.time() - start)/60:.2f} minutes")

        print("\n--- Extracting hematoxylin ---")
        start = time.time()
        moving_arr = stains.hematoxylin_deconvolution(moving_arr)
        # The array is 2D now, so the axes string and the axis indices to be recomputed
        axes = "YX"
        y_axis_idx, x_axis_idx = io_utils.spatial_axes(axes)
        moving_xml = io_utils.single_channel_metadata(moving_xml, ome_type=moving_arr.dtype.name, name="Hematoxylin")
        print(f"shape {moving_arr.shape}, dtype {moving_arr.dtype}, axes '{axes}' (Y={y_axis_idx}, X={x_axis_idx})")
        print(f"extracted in {(time.time() - start)/60:.2f} minutes")

    # Stage the transformed image (LZW-compressed OME-TIFF)
    print("\n--- Staging images (for VALIS to use) ---")
    print(f"staged dir: {staged_dir}")
    print(f"moving    : {moving_name}")
    print(f"reference : {reference_name}")
    start = time.time()
    valis_runner.stage_images(
        moving_arr=moving_arr,
        moving_xml=moving_xml,
        moving_name=moving_name,
        reference_name=reference_name,
        reference_path=fixed_image,
        staged_dir=staged_dir,
        compression=compression,
    )
    print(f"staged in {(time.time() - start)/60:.2f} minutes")

    # Run VALIS registration on the staged images. 
    # Hold the JVM open when there is still an RGB image to warp (H&E case)
    print("\n--- Running registration ---")
    registered_path = valis_runner.run_registration(
        moving=staged_dir / moving_name,
        fixed=staged_dir / reference_name,
        output_dir=output_dir,
        kill_jvm=rgb_staged_path is None,
    )

    # Apply the transforms VALIS just learned on the hematoxylin to the RGB (original) H&E image
    rgb_registered_path = None
    if rgb_staged_path is not None:
        print("\n--- Warping the RGB with the same transforms ---")
        start = time.time()
        pickle_path = valis_runner.find_registrar_pickle(output_dir / "valis_results")
        print(f"registrar: {pickle_path}")
        rgb_registered_path = output_dir / "registered_slides" / f"{rgb_staged_path.name.split('.')[0]}_registered.ome.tiff"
        valis_runner.register_pickle(
            was_registered=staged_dir / moving_name,
            to_register=rgb_staged_path,
            pickle_path=pickle_path,
            output_path=rgb_registered_path,
        )
        print(f"warped in {(time.time() - start)/60:.2f} minutes")

    print("\n--- Done ---")
    print(f"registered slide: {registered_path}")
    if rgb_registered_path is not None:
        print(f"registered RGB  : {rgb_registered_path}")
    print(f"outputs in {output_dir}")


if __name__ == "__main__":
    main()