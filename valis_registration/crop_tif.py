"""
Crops an OME-TIFF to a specified region defined
by coordinates (x, y, width, height).

Preserves OME-TIFF metadata and structure.
The cropped OME-TIFF is LZW compressed and 1024x1024 tiled.

Inputs:
- input_image: Path to the input OME-TIFF
- output_image: Path to write the cropped OME-TIFF
- x: X-coordinate of the top-left corner of the crop region
- y: Y-coordinate of the top-left corner of the crop region
- width: Width of the crop region (px)
- height: Height of the crop region (px)

Submit with 
python crop_tif.py --input_image path/to/input.ome.tiff --output_image path/to/output.ome.tiff --x 0 --y 0 --width 1024 --height 1024
"""

import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

import io_utils

pyramid_ns = "openmicroscopy.org/PyramidResolution"

def crop_metadata(xml, width, height):
    """
    Rewrite the source OME-XML to describe the cropped output.
    Returns the modified OME-XML as a string.
    """
    if xml is None:
        return None
    root = ET.fromstring(xml)
    ns = {"ome": io_utils.ome_ns}

    # SizeX/SizeY become the crop dimensions
    for pixels in root.findall(".//ome:Pixels", ns):
        pixels.set("SizeX", str(width))
        pixels.set("SizeY", str(height))

    # Remove the pyramid annotation, since the cropped image is a single resolution level
    for annotations in root.findall("ome:StructuredAnnotations", ns):
        for annotation in annotations.findall("ome:MapAnnotation", ns):
            if annotation.get("Namespace") == pyramid_ns:
                annotations.remove(annotation)
        if len(annotations) == 0:
            root.remove(annotations)

    return io_utils.ome_tostring(root)

def parse_args():
    parser = argparse.ArgumentParser(description="Crop an OME-TIFF to a region")
    parser.add_argument("--input_image", type=Path, required=True, help="Path to the input OME-TIFF")
    parser.add_argument("--output_image", type=Path, required=True, help="Path to write the cropped OME-TIFF")
    parser.add_argument("--x", type=int, required=True, help="X-coordinate of the top-left corner of the crop region")
    parser.add_argument("--y", type=int, required=True, help="Y-coordinate of the top-left corner of the crop region")
    parser.add_argument("--width", type=int, required=True, help="Width of the crop region")
    parser.add_argument("--height", type=int, required=True, help="Height of the crop region")
    return parser.parse_args()

def main():
    args = parse_args()

    print("\n--- Inputs ---")
    print(f"input image : {args.input_image}")
    print(f"output image: {args.output_image}")
    print(f"crop region : x={args.x}, y={args.y}, width={args.width}, height={args.height}")

    arr, xml, axes = io_utils.read_ome(args.input_image)
    y_axis_idx, x_axis_idx = io_utils.spatial_axes(axes)
    print(f"shape {arr.shape}, dtype {arr.dtype}, axes '{axes}' (Y={y_axis_idx}, X={x_axis_idx})")

    # Build the slicing object to crop the array
    crop = [slice(None)] * arr.ndim
    crop[y_axis_idx] = slice(args.y, args.y + args.height)
    crop[x_axis_idx] = slice(args.x, args.x + args.width)
    cropped_arr = arr[tuple(crop)]
    print(f"cropped shape {cropped_arr.shape}")

    # Slicing clamps at the image edge, so the metadata takes the actual shape
    cropped_xml = crop_metadata(xml, cropped_arr.shape[x_axis_idx], cropped_arr.shape[y_axis_idx])

    print("\n--- Writing cropped image ---")
    io_utils.write_ome_tif(str(args.output_image), cropped_arr, cropped_xml)
    print(f"wrote {args.output_image}")

if __name__ == "__main__":
    main()
