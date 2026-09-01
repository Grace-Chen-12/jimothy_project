"""
Goal: Convert nuclei label mask(s) (Zarr) into a Geojson annotation file that can be imported directly into QuPath.

Inputs:
- Nuclei segmentation mask Zarr file (required)
- Micronuclei segmentation mask Zarr file (optional)
Output:
- GeoJSON file with QuPath annotations, one Polygon per nucleus, classified as "Nucleus" or "Micronucleus" based on the optional micronuclei mask.
"""

import numpy as np
import zarr
import cv2
import shapely
from shapely.geometry import Polygon, mapping
from shapely.affinity import translate as shp_translate
from shapely.ops import unary_union
from skimage.measure import regionprops
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--nuclei_mask', type=str, required=True,
                     help='Path to the full nuclei label mask (.zarr) from classify_all_nuclei.py --output')
parser.add_argument('--micronuclei_mask', type=str, default=None,
                     help='Path to the micronuclei-only label mask (.zarr) from filter_for_mn.py')
parser.add_argument('--output', type=str, default=None, help='Path to write the GeoJSON annotations')
parser.add_argument('--tile_size', type=int, default=4096, help='Must match the value used in classify_all_nuclei.py')
parser.add_argument('--overlap', type=int, default=256, help='Must match the value used in classify_all_nuclei.py')
parser.add_argument('--simplify', type=float, default=0.0,
                     help='cv2.approxPolyDP epsilon in pixels (0 disables simplification, the default). ')
args = parser.parse_args()

NUCLEI_MASK = args.nuclei_mask
MN_MASK = args.micronuclei_mask
OUTPUT = args.output or str(Path(NUCLEI_MASK).with_suffix('')) + '.geojson'
TILE_SIZE = args.tile_size
OVERLAP = args.overlap
SIMPLIFY = args.simplify

CLASS_COLORS = {
    'Nucleus': [31, 119, 180],
    'Micronucleus': [214, 39, 40],
}


def make_tiles(height, width, tile_size, overlap):
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


def _polygonal_parts(geom):
    """Keep only the polygonal (nonzero-area) components of geom, discarding spray points. 
    Returns a Polygon/MultiPolygon, or None if nothing is left."""
    if geom.is_empty:
        return None
    if geom.geom_type in ('Polygon', 'MultiPolygon'):
        return geom
    if geom.geom_type == 'GeometryCollection':
        polys = [g for g in geom.geoms if g.geom_type in ('Polygon', 'MultiPolygon') and not g.is_empty]
        return unary_union(polys) if polys else None
    return None


def mask_to_geometry(local_mask, simplify_eps):
    """local_mask: 2D bool array (prop.image). Returns a shapely
    Polygon/MultiPolygon in local_mask coordinates, or None if no valid
    geometry was found."""
    contours, _ = cv2.findContours(
        local_mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
    )
    if not contours:
        return None

    if simplify_eps > 0:
        contours = [cv2.approxPolyDP(c, simplify_eps, closed=True) for c in contours]

    polys = []
    for c in contours:
        pts = c.reshape(-1, 2)
        if len(pts) < 3:
            continue
        poly = Polygon(pts)
        if not poly.is_valid:
            poly = _polygonal_parts(shapely.make_valid(poly))
        if poly is not None and not poly.is_empty:
            polys.append(poly)

    if not polys:
        return None
    geom = polys[0] if len(polys) == 1 else unary_union(polys)
    return None if geom.is_empty else geom


nuclei_store = zarr.open(NUCLEI_MASK, mode='r')
HEIGHT, WIDTH = nuclei_store.shape
print(f"Mask shape (Y, X): {nuclei_store.shape}")

mn_store = None
if MN_MASK:
    mn_store = zarr.open(MN_MASK, mode='r')
    assert mn_store.shape == (HEIGHT, WIDTH), (
        f"Micronuclei mask shape {mn_store.shape} != nuclei mask shape {(HEIGHT, WIDTH)}"
    )

tiles = list(make_tiles(HEIGHT, WIDTH, TILE_SIZE, OVERLAP))
print(f"Processing {len(tiles)} tiles of size {TILE_SIZE}px with {OVERLAP}px overlap")

class_counts = {}
n_features = 0

with open(OUTPUT, 'w') as out:
    out.write('{"type": "FeatureCollection", "features": [\n')
    first_feature = True

    for i, (y0, y1, x0, x1) in enumerate(tiles):
        nuc_tile = np.asarray(nuclei_store[y0:y1, x0:x1])
        if not nuc_tile.any():
            continue
        mn_tile = np.asarray(mn_store[y0:y1, x0:x1]) if mn_store is not None else None

        core_y0 = OVERLAP // 2 if y0 > 0 else 0
        core_x0 = OVERLAP // 2 if x0 > 0 else 0
        core_y1 = (y1 - y0) - (OVERLAP // 2 if y1 < HEIGHT else 0)
        core_x1 = (x1 - x0) - (OVERLAP // 2 if x1 < WIDTH else 0)

        n_kept = 0
        for prop in regionprops(nuc_tile):
            cy, cx = prop.centroid
            if not (core_y0 <= cy < core_y1 and core_x0 <= cx < core_x1):
                continue  # centroid belongs to a neighboring tile

            if mn_tile is not None:
                rows, cols = prop.coords[:, 0], prop.coords[:, 1]
                cls_name = 'Micronucleus' if mn_tile[rows, cols].any() else 'Nucleus'
            else:
                cls_name = 'Nucleus'

            row_min, col_min = prop.bbox[0], prop.bbox[1]
            geom = mask_to_geometry(prop.image, SIMPLIFY)
            if geom is None:
                continue
            geom = shp_translate(geom, xoff=col_min + x0, yoff=row_min + y0)

            feature = {
                "type": "Feature",
                "geometry": mapping(geom),
                "properties": {
                    "objectType": "annotation",
                    "classification": {"name": cls_name, "color": CLASS_COLORS[cls_name]},
                    "name": f"Label {prop.label}",
                },
            }

            if not first_feature:
                out.write(',\n')
            out.write(json.dumps(feature))
            first_feature = False

            n_kept += 1
            n_features += 1
            class_counts[cls_name] = class_counts.get(cls_name, 0) + 1

        print(f"Tile {i + 1}/{len(tiles)} ({y0}:{y1}, {x0}:{x1}): {n_kept} objects written")

    out.write('\n]}\n')

print(f"\n{n_features} annotations written to {OUTPUT}")
for cls_name, count in sorted(class_counts.items()):
    print(f"  {cls_name}: {count}")
