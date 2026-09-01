## VALIS Registration

Apply [VALIS](https://www.nature.com/articles/s41467-023-40218-9) to help register ome.tiff whole-slide images (WSI) on BioHPC. These scripts live in `valis_registration/`.

### Notes

The moving image may be an H&E or IF image. The reference image must be an IF image. For all IF images, it is expected they contain a DAPI channel and that the DAPI channel is the first channel (channel 0).

### Installation/Set-up

Note that BioHPC uses `docker1` instead of `docker`.

Pull the docker image:

```
docker1 pull cdgatenbee/valis-wsi:1.2.0
```

Check the name of the pulled docker image:

```
docker1 images
```

All outputs written from inside a container are owned by root on the host. So, reclaim access to VALIS outputs with

```
docker1 claim /workdir/$USER/whatever_dir/
```

### Cropping Images

Sometimes, you want to crop images before running registration.

It may be easier to crop ome.tiff images with QuPath. Use the rectangle shape tool to select the section of the image you would like to crop. Then, select _File > Export Images > OME-TIFF_. I recommend the following settings (the QuPath default):

- Compression type: LZW (loseless)
- Pyramidal scale: 4.0 (each pyramid image is downsampled by a factor of 4)
- Tile size: 256 px
- Parallelize export: Keep setting checked

Then, press "OK" and wait for the file export process to complete.

If you don't have QuPath, you can use `crop_tif.py` to crop images.
Run with

```
python crop_tif.py --input_image path/to/input.ome.tiff --output_image path/to/output.ome.tiff --x 0 --y 0 --width 1024 --height 1024
```

in which x is the X-coordinate of the top-left corner of the crop region, y is the Y-coordinate of the top-left corner of the crop region, width is the width of the crop region in px, and height is the height of the crop region in px.

### Running Registration

Run registration with _register.sh_. An example job submission is

```
sbatch register.sh --moving_image "/workdir/$USER/path/to/moving.ome.tif" --fixed_image "/workdir/$USER/path/to/fixed.ome.tif" --output_dir /workdir/$USER/output/dir --moving_modality HE --transform rot90
```

Change the --mail-user in the sbatch directives to your own email before running _register.sh_.

The options for moving_modality are "HE" (H&E) and "IF". Note that the reference image is always IF. The main difference between the two options is that channel axes are handled differently, and a hematoxylin deconvolution is performed if H&E is selected.

Use the --transform tag to apply a transform to the moving image before registration (VALIS cannot perform reflections). The options for --transform are identity, rot90, rot180, rot270, flipud, fliplr, transpose, and transverse. All rotations are counterclockwise. "transpose" describes a fliplr then rot90. "transverse" describes a flipud then rot90.

The output directory will contain the following directories: `registered_slides`, `valis_results`, and `staged`. `registered_slides` contains the slide registered with VALIS. `valis_results` contains various files, including a summary spreadsheet of alignment results, a pickled version of the registrar, and image thumbnails that can be used for quick inspection. `staged` contains the two processed ome.tiff images on which registration was actually performed.

When _--moving_modality_ is HE, the output directory will also contain a `staged_rgb/` directory and a second registered file in `registered_slides` (the RGB version of the registered slide, in addition to the deconvolved hematoxylin version). The transforms VALIS learns on the hematoxylin channel are re-applied to the RGB image.

### Xenium-Explorer Compatible Output

The VALIS outputs are not directly compatible with Xenium Explorer. Use _make_pyramid.sh_ to make a Xenium-compatible version of a tiff. Change the --mail-user in the sbatch directives to your own email before running _make_pyramid.sh_.

An example job submission for _make_pyramid.sh_ is

```
sbatch make_pyramid.sh /workdir/$USER/path/to/image.ome.tiff
```

The output is written next to the input file, named ..._xenium_compat.ome.tif_.

### Deleting Valis Outputs

You will need to reclaim VALIS outputs or you will be unable to delete them. Reclaim with

```
docker1 claim /workdir/$USER/whatever_dir/
```

## Micronuclei Segmentation

These scripts live in `mn_segmentation/`. Before running them, change the `--mail-user` in the sbatch directives to your own email and adjust the `source .../miniconda3/etc/profile.d/conda.sh` line to match your conda install path. _classify_all_nuclei.sh_ segments all nuclei (PN and MN). _filter_for_mn.sh_ then classifies the segmented nuclei into primary nuclei (PN) and micronuclei (MN). Segmentation masks are provided in a .zarr format. The output .zarr files can be turned into annotations that can be imported into QuPath (.geojson format) with _zarr_to_qupath_annotations.sh_.

### Install `cellpose` with `pip`

If you haven't already, setup Conda on BioHPC. This is a helpful resource for installation: https://biohpc.cornell.edu/lab/doc/Software_Installation_exercises1.html

Create the env and install:

```
conda create -n cellpose python=3.11 -y
conda activate cellpose
pip install cellpose
```

### Segment All Nuclei

Submit the _classify_all_nuclei.sh_ script to segment all nuclei (PN and MN). Submit with
`sbatch classify_all_nuclei.sh /path/to/ch0000_dapi.ome.tif /path/to/output_mask.zarr channel_min channel_max model` in which

- _/path/to/ch0000_dapi.ome.tif_ is the file path to the _ch0000_dapi.ome.tif_ file from the morphology_focus folder from the Xenium outputs (required)
- _/path/to/output_mask.zarr_ is the path to write the .zarr nuclei masks (optional, will be written in same DIR as ch0000_dapi.ome.tif file if not provided)
- `channel_min` is a float representing DAPI normalization floor (optional, suitable value will be estimated automatically if not provided)
- `channel_max` is a float representing DAPI normalization ceiling (optional, suitable value will be estimated automatically if not provided)
- `channel_min` and `channel_max` must be given together. If only one is provided it is ignored and both are estimated from the image.
- `model` is the cellpose pre-trained model to use (optional, string, default is `cpsam_v2`). You can provide a path to your own pre-trained model if you have one.

Arguments are positional and must be given in the order above (the scripts add the `--` tags themselves). To skip a middle argument, you have to pass an empty string "" as a placeholder.

### Filter For Micronuclei

Submit the _filter_for_mn.sh_ script to classify which of the nuclei found in the 1st step are MN. Submit with `sbatch filter_for_mn.sh /path/to/nuclei_mask.zarr /path/to/cells.zarr.zip /path/to/output /path/to/post_if.zarr` in which:

- _/path/to/nuclei_mask.zarr_ is the output from the 1st step (required)
- _/path/to/cells.zarr.zip_ is the _cells.zarr.zip_ file from the Xenium outputs (required)
- _/path/to/output_ is the base path the outputs are named from (optional, will be written in same DIR as nuclei_mask if not provided)
- _/path/to/post_if.zarr_ is the registered post-IF .zarr cache, shape (n_channels, Y, X) (optional, per-nucleus cGAS intensity columns are added to the summary CSV when provided)

Three files are written from that base path: _<base>\_micronuclei_mask.zarr_ (the MN-only label mask), _<base>\_micronuclei_summary.csv_ (per-nucleus area, circularity, centroid, cell assignment, classification and exclusion reason), and _<base>\_all_nuclei.csv_ (nucleus_label, cell_id, is_single, is_micronucleus, area_px, cell_boundary_method).

### Convert Zarr to QuPath Annotations

Submit the _zarr_to_qupath_annotations.sh_ script to convert the Zarr files to a geojson file that can be read by QuPath. Submit with `sbatch zarr_to_qupath_annotations.sh /path/to/nuclei_mask.zarr /path/to/micronuclei_mask.zarr /path/to/output.geojson` in which

- _/path/to/nuclei_mask.zarr_ is the output from the 1st step (required)
- _/path/to/micronuclei_mask.zarr_ is the output from the 2nd step (optional, only nuclei annotations will be written if this is not provided)
- _/path/to/output.geojson_ is the path in which output geojson is written (optional, same DIR as nuclei_mask if not provided)

## Dependencies

- [tifffile](https://github.com/cgohlke/tifffile)
- [numpy](https://numpy.org/)
- [opencv](https://opencv.org/)
- [skimage](https://scikit-image.org/)
- [scipy](https://scipy.org/)
- [zarr](https://zarr.readthedocs.io/)
- [shapely](https://shapely.readthedocs.io/)
- [cellpose](https://github.com/MouseLand/cellpose)
