"""
Goal: classify which nuclei (from existing nuclei mask produced by classify_all_nuclei.py) are MN.
A nucleus is classified as a MN if it shares a Xenium cell boundary with a larger nucleus (the PN) in the same cell.

The following filters were added to reduce spurious MN annotations:
- The candidate MN must be within 1/256-3/5 of the PN area. (A nucleus outside that size ratio range is considered neither PN nor a MN, and is left out of the MN mask.)
- The candidate MN must not touch any primary nucleus or any larger nucleus in the nuclei mask.
- The candidate MN must be within a cell boundary determined by the 'Segmented by nucleus expansion of 5.0um' segmentation method.

INPUTS:
- Zarr file with nuclei mask (produced by classify_all_nuclei.py)
- cells.zarr.zip Xenium file

OUTPUTS:
- Zarr file with MN mask
- CSV file with per-nucleus summary (including area, circularity, centroid and, if --post_if_zarr is given, mean/median/max cGAS intensity) for downstream use
- CSV file with all-nuclei summary (nucleus_label, cell_id, is_single, is_micronucleus, area_px, cell_boundary_method)
"""

import argparse
import csv
from pathlib import Path
import numpy as np
import zarr
from scipy import ndimage as ndi
from skimage.measure import perimeter_crofton

def open_zarr(path: str) -> zarr.Group:
    '''Open a Zarr file from path'''
    store = (
        zarr.storage.ZipStore(path, mode="r") if path.endswith(".zip")
        else zarr.storage.DirectoryStore(path)
    )
    return zarr.open_group(store=store, mode="r")


parser = argparse.ArgumentParser()
parser.add_argument('--nuclei_mask', type=str, required=True, help='Path to the nuclei .zarr mask produced by classify_all_nuclei.py')
parser.add_argument('--cells_zarr', type=str, required=True, help='Path to the Xenium cells.zarr.zip file')
parser.add_argument('--output', type=str, default=None, help='Base path for outputs (defaults to next to --nuclei_mask)')
# registered post-IF cache, shape (n_channels, Y, X), same (Y, X) extent as DAPI/masks above
parser.add_argument('--post_if_zarr', type=str, required=False, help='Path to the registered post-IF .zarr file')
args = parser.parse_args()

NUCLEI_MASK_PATH = args.nuclei_mask
CELLS_ZARR = args.cells_zarr
OUTPUT_BASE = args.output or str(Path(NUCLEI_MASK_PATH).with_suffix(''))
POST_IF_ZARR = args.post_if_zarr

# Cell-boundary segmentation methods matching this substring 
# (case-insensitive) are excluded from MN classification
REJECTED_METHOD_SUBSTR = ['nucleus expansion']

# Broad ratio bounds for MN as a basic QC filter
MN_MIN_SIZE_RATIO = 1 / 256
MN_MAX_SIZE_RATIO = 3 / 5

MN_PN_CONTACT_DILATION = 1
CONTACT_STRUCTURE = np.ones((3, 3), dtype=bool)


mask_store = zarr.open(NUCLEI_MASK_PATH, mode='r')
HEIGHT, WIDTH = mask_store.shape
print(f"Nuclei mask shape: {mask_store.shape}")

cells_zip_root = open_zarr(CELLS_ZARR)
print("Cells zip root info:", cells_zip_root.info)
print("Cells zip root structure:", cells_zip_root.tree())

# mask_index=1 is the cell segmentation mask
cellseg_zarr = cells_zip_root["masks"][1]

print("Cell segmentation mask shape (should match the nuclei mask):", cellseg_zarr.shape)
print("Total cells detected (cell_id count):", cells_zip_root["cell_id"].shape[0])

assert cellseg_zarr.shape == (HEIGHT, WIDTH), (
    f"Xenium cell mask shape {cellseg_zarr.shape} != nuclei mask shape {(HEIGHT, WIDTH)}. "
    "The morphology image and cells.zarr.zip must share the same pixel grid "
)

segmentation_methods = list(cells_zip_root.attrs["segmentation_methods"])
rejected_method_ids = {
    i for i, name in enumerate(segmentation_methods)
    if any(substr in name.lower() for substr in REJECTED_METHOD_SUBSTR)
}
print(f"Rejecting cells whose boundary used: {[segmentation_methods[i] for i in rejected_method_ids]}")

# polygon_sets_index=1 contains cell segmentation polygons
cell_polygons = cells_zip_root["polygon_sets"][1]
poly_cell_index = np.asarray(cell_polygons["cell_index"])
poly_method = np.asarray(cell_polygons["method"])
method_by_cell_label = np.full(int(poly_cell_index.max()) + 1, -1, dtype=np.int32)
method_by_cell_label[poly_cell_index] = poly_method

def cell_method_name(cell_id: int) -> str:
    '''Name of the segmentation method used for a given Xenium cell label, or '' if unknown'''
    if cell_id <= 0 or cell_id >= len(method_by_cell_label):
        return ''
    method_id = method_by_cell_label[cell_id]
    return segmentation_methods[method_id] if method_id >= 0 else ''

if POST_IF_ZARR:
    CGAS_CHANNEL = 1  # C2 of post-IF cache
    post_if = zarr.open(POST_IF_ZARR, mode='r')
    print("Registered post-IF cache shape (n_channels, Y, X):", post_if.shape)
    assert post_if.shape[1:] == (HEIGHT, WIDTH), (
        f"Registered post-IF cache shape {post_if.shape} != nuclei mask shape {(HEIGHT, WIDTH)}. "
        "The morphology image and registered post-IF cache must share the same pixel grid "
    )

# Loads the full nuclei label mask into RAM
print("Loading full nuclei mask for classification...")
nuclei_full = np.asarray(mask_store)


print("Computing per-nucleus properties (area, centroid)...")

# Obtain row/column indices for all non-zero pixels in nuclei mask
fg_rows, fg_cols = np.nonzero(nuclei_full)
# Obtain which nucleus (label) each pixel belongs to
fg_labels = nuclei_full[fg_rows, fg_cols]

max_label = int(fg_labels.max())
area_by_label = np.bincount(fg_labels, minlength=max_label + 1)
row_sum_by_label = np.bincount(fg_labels, weights=fg_rows.astype(np.float64), minlength=max_label + 1)
col_sum_by_label = np.bincount(fg_labels, weights=fg_cols.astype(np.float64), minlength=max_label + 1)

# Summarise the cGAS channel over each nucleus' pixels while the foreground index arrays are still around
nucleus_mean_cgas = None  # label -> mean cGAS intensity
cgas_median = None        # label -> median cGAS intensity
cgas_max = None           # label -> max cGAS intensity
if POST_IF_ZARR:
    print("Summarising cGAS intensity per nucleus...")
    cgas_plane = np.asarray(post_if[CGAS_CHANNEL])
    cgas_values = cgas_plane[fg_rows, fg_cols].astype(np.float64)
    del cgas_plane
    cgas_sum_by_label = np.bincount(fg_labels, weights=cgas_values, minlength=max_label + 1)

    # Sort pixels by (label, intensity) so each nucleus' values sit in one sorted run,
    # which makes the per-nucleus median and max plain index lookups.
    order = np.lexsort((cgas_values, fg_labels))
    sorted_labels = fg_labels[order]
    sorted_values = cgas_values[order]
    del order, cgas_values
    run_labels, run_starts, run_counts = np.unique(
        sorted_labels, return_index=True, return_counts=True
    )
    mid = run_starts + run_counts // 2
    median_values = np.where(
        run_counts % 2 == 1,
        sorted_values[mid],
        (sorted_values[mid] + sorted_values[np.maximum(mid - 1, run_starts)]) / 2,
    )
    cgas_median = dict(zip(run_labels.tolist(), median_values.tolist()))
    cgas_max = dict(zip(run_labels.tolist(), sorted_values[run_starts + run_counts - 1].tolist()))
    del sorted_labels, sorted_values

    nucleus_mean_cgas = dict(zip(
        run_labels.tolist(), (cgas_sum_by_label[run_labels] / area_by_label[run_labels]).tolist()
    ))
    del cgas_sum_by_label

del fg_rows, fg_cols, fg_labels

labels = np.nonzero(area_by_label)[0]
labels = labels[labels != 0]
areas = area_by_label[labels]
centroid_row = row_sum_by_label[labels] / areas
centroid_col = col_sum_by_label[labels] / areas
print(f"  {len(labels)} nuclei total")


print("Computing per-nucleus circularity...")

# Circularity = 4*pi*area / perimeter^2; 1.0 for a perfect disk; -> 0 for elongated
# Perimeter uses skimage's Crofton estimator

object_slices = ndi.find_objects(nuclei_full)
perimeter_by_label = np.zeros(len(labels), dtype=np.float64)
for i, lbl in enumerate(labels.tolist()):
    bbox = object_slices[lbl - 1]
    if bbox is None:
        continue
    # Pad by 1 so the nucleus never touches the crop edge, which would clip its boundary
    perimeter_by_label[i] = perimeter_crofton(np.pad(nuclei_full[bbox] == lbl, 1))

with np.errstate(divide='ignore', invalid='ignore'):
    circularity = np.where(
        perimeter_by_label > 0, 4 * np.pi * areas / perimeter_by_label**2, np.nan
    )
nucleus_circularity = dict(zip(labels.tolist(), circularity.tolist()))
n_degenerate = int((~np.isfinite(circularity)).sum())
if n_degenerate:
    print(f"  {n_degenerate} nuclei have no measurable perimeter; circularity left blank")

def circularity_str(lbl: int) -> str:
    '''Circularity of a nucleus rounded for CSV output, or '' if undefined'''
    value = nucleus_circularity[lbl]
    return '' if not np.isfinite(value) else f'{value:.5f}'

def cgas_cells(lbl: int) -> list:
    '''Trailing CSV cells holding the nucleus' mean/median/max cGAS intensity, empty when no post-IF cache was given'''
    if nucleus_mean_cgas is None:
        return []
    return [
        f'{nucleus_mean_cgas[lbl]:.4f}',
        f'{cgas_median[lbl]:.4f}',
        f'{cgas_max[lbl]:.4f}',
    ]

nucleus_centroid = dict(zip(
    labels.tolist(), zip(centroid_row.tolist(), centroid_col.tolist())
))

# Assign each nucleus to Xenium cell containing its centroid pixel
# Nuclei whose centroid lands outside all cells are left out of classification entirely
row_idx = np.round(centroid_row).astype(np.int64)
col_idx = np.round(centroid_col).astype(np.int64)
cell_ids = cellseg_zarr.vindex[row_idx, col_idx]

nucleus_area = dict(zip(labels.tolist(), areas.tolist()))

cell_to_nuclei = {}

all_cell_to_nuclei = {}
excluded_labels = {}  # label -> exclusion reason
rejected_cell_count = 0
for lbl, cell_id in zip(labels.tolist(), cell_ids.tolist()):
    if cell_id <= 0:
        excluded_labels[lbl] = 'no_cell_match'
        continue
    all_cell_to_nuclei.setdefault(cell_id, []).append(lbl)
    method_id = method_by_cell_label[cell_id] if cell_id < len(method_by_cell_label) else -1
    if method_id in rejected_method_ids:
        excluded_labels[lbl] = 'nucleus_expansion_boundary'
        rejected_cell_count += 1
        continue
    cell_to_nuclei.setdefault(cell_id, []).append(lbl)

print(f"{len(excluded_labels)} nuclei excluded from classification "
      f"(of {len(labels)} total nuclei, {rejected_cell_count} cells rejected for nucleus expansion boundary segmentation method)")

# The largest nucleus in each cell, over every cell rather than only the eligible ones
primary_by_cell = {
    cell_id: max(nuclei_in_cell, key=lambda lbl: nucleus_area[lbl])
    for cell_id, nuclei_in_cell in all_cell_to_nuclei.items()
}
primary_labels = set(primary_by_cell.values())
print(f"{len(primary_labels)} nuclei are the primary of their cell")


def touching_labels(candidate_label: int) -> list:
    '''
    Every other nucleus within MN_PN_CONTACT_DILATION px of candidate_label.
    '''
    bbox = object_slices[candidate_label - 1]
    if bbox is None:
        return []
    pad = MN_PN_CONTACT_DILATION
    rows = slice(max(bbox[0].start - pad, 0), min(bbox[0].stop + pad, HEIGHT))
    cols = slice(max(bbox[1].start - pad, 0), min(bbox[1].stop + pad, WIDTH))
    crop = nuclei_full[rows, cols]
    reach = ndi.binary_dilation(
        crop == candidate_label, structure=CONTACT_STRUCTURE, iterations=pad
    )
    # reach covers the candidate as well as its 1 px collar, so drop background and the
    # candidate's own label from what the collar picks up.
    return [
        int(other) for other in np.unique(crop[reach])
        if other != 0 and other != candidate_label
    ]


micronucleus_labels = []
classification_by_label = {}  # label -> (cell_id, primary_label, primary_area, is_primary, is_micronucleus, reason)
n_size_rejected = 0
n_primary_contact_rejected = 0
n_larger_contact_rejected = 0
for cell_id, nuclei_in_cell in cell_to_nuclei.items():
    primary_label = primary_by_cell[cell_id]
    primary_area = nucleus_area[primary_label]
    for lbl in nuclei_in_cell:
        is_primary = lbl == primary_label
        if is_primary:
            classification_by_label[lbl] = (cell_id, primary_label, primary_area, True, False, '')
            continue
        size_ratio = nucleus_area[lbl] / primary_area
        is_micronucleus = MN_MIN_SIZE_RATIO <= size_ratio <= MN_MAX_SIZE_RATIO
        if not is_micronucleus:
            n_size_rejected += 1
            reason = 'mn_size_ratio_out_of_range'
        else:
            touching = touching_labels(lbl)

            if any(other in primary_labels for other in touching):
                is_micronucleus = False
                n_primary_contact_rejected += 1
                reason = 'touches_primary_nucleus'

            elif any(nucleus_area[other] > nucleus_area[lbl] for other in touching):
                is_micronucleus = False
                n_larger_contact_rejected += 1
                reason = 'touches_larger_nucleus'
            else:
                micronucleus_labels.append(lbl)
                reason = ''
        classification_by_label[lbl] = (cell_id, primary_label, primary_area, False, is_micronucleus, reason)

print(f"{len(micronucleus_labels)} nuclei classified as micronuclei "
      f"(of {len(nucleus_area)} total nuclei, {len(cell_to_nuclei)} "
      f"cells matched to >=1 eligible nucleus, {n_size_rejected} non-primary nuclei rejected "
      f"for falling outside the [{MN_MIN_SIZE_RATIO:.4f}, {MN_MAX_SIZE_RATIO:.4f}] PN-size ratio, "
      f"{n_primary_contact_rejected} rejected for touching a primary nucleus, "
      f"{n_larger_contact_rejected} for touching a larger non-primary nucleus)")

SUMMARY_PATH = OUTPUT_BASE + '_micronuclei_summary.csv'
with open(SUMMARY_PATH, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow([
        'nucleus_label', 'area_px', 'circularity', 'centroid_row', 'centroid_col',
        'cell_id', 'cell_segmentation_method',
        'primary_label', 'primary_area_px', 'is_primary', 'is_micronucleus',
        'exclusion_reason',
    ] + ([] if nucleus_mean_cgas is None else
         ['mean_cgas_intensity', 'median_cgas_intensity', 'max_cgas_intensity']))
    for lbl, cell_id in zip(labels.tolist(), cell_ids.tolist()):
        area = nucleus_area[lbl]
        c_row, c_col = nucleus_centroid[lbl]
        if lbl in classification_by_label:
            c_id, primary_label, primary_area, is_primary, is_micronucleus, reason = classification_by_label[lbl]
            writer.writerow([
                lbl, area, circularity_str(lbl), f'{c_row:.2f}', f'{c_col:.2f}',
                c_id, cell_method_name(c_id),
                primary_label, primary_area, is_primary, is_micronucleus, reason,
            ] + cgas_cells(lbl))
        else:
            writer.writerow([
                lbl, area, circularity_str(lbl), f'{c_row:.2f}', f'{c_col:.2f}',
                cell_id if cell_id > 0 else '', cell_method_name(cell_id),
                '', '', False, False, excluded_labels[lbl],
            ] + cgas_cells(lbl))
print(f"Micronucleus summary saved to {SUMMARY_PATH}")

ALL_NUCLEI_PATH = OUTPUT_BASE + '_all_nuclei.csv'
with open(ALL_NUCLEI_PATH, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow([
        'nucleus_label', 'cell_id', 'is_single', 'is_micronucleus', 'area_px', 'cell_boundary_method',
    ])
    for lbl, cell_id in zip(labels.tolist(), cell_ids.tolist()):
        is_micronucleus = classification_by_label.get(lbl, (None, None, None, None, False, ''))[4]
        if cell_id > 0:
            is_single = len(all_cell_to_nuclei[cell_id]) == 1
        else:
            is_single = ''
        writer.writerow([
            lbl, cell_id if cell_id > 0 else '', is_single, is_micronucleus,
            nucleus_area[lbl], cell_method_name(cell_id),
        ])
print(f"All-nuclei summary saved to {ALL_NUCLEI_PATH}")

MN_MASK_PATH = OUTPUT_BASE + '_micronuclei_mask.zarr'
mn_mask = zarr.open(
    MN_MASK_PATH, mode='w', shape=(HEIGHT, WIDTH),
    chunks=mask_store.chunks, dtype=np.uint32,
)

# Zero out everything except the classified micronuclei in place
nuclei_full[~np.isin(nuclei_full, micronucleus_labels)] = 0
mn_mask[:] = nuclei_full
print(f"Micronuclei mask saved to {MN_MASK_PATH}")
