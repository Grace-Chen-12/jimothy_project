import tifffile
import xml.etree.ElementTree as ET

ome_ns = "http://www.openmicroscopy.org/Schemas/OME/2016-06"
xy_pairs = (("SizeX", "SizeY"), ("PhysicalSizeX", "PhysicalSizeY"), ("PhysicalSizeXUnit", "PhysicalSizeYUnit"))

micron = "\u00b5m"
cm_per_unit = {
    "m": 100.0, "cm": 1.0, "mm": 0.1,
    micron: 1e-4, "\u03bcm": 1e-4, "um": 1e-4,
    "nm": 1e-7, "pm": 1e-10, "\u00c5": 1e-8,
}

def resolution_from_xml(xml):
    """
    Return image resolution as (res_x, res_y) in px/cm, or None if unavailable.
    """
    if xml is None:
        return None
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return None

    # all dimensional data is under the Pixels element
    pixels = root.find(f".//{{{ome_ns}}}Pixels")
    if pixels is None:
        return None
    
    resolution = []
    for size_attr, unit_attr in (("PhysicalSizeX", "PhysicalSizeXUnit"),
                                 ("PhysicalSizeY", "PhysicalSizeYUnit")):
        # Note the OME default when the unit is absent is micrometres
        cm = cm_per_unit.get(pixels.get(unit_attr, micron))
        try:
            size = float(pixels.get(size_attr))
        except (TypeError, ValueError):
            return None
        if cm is None or size <= 0:
            return None
        resolution.append(1.0 / (size * cm))
    return tuple(resolution)

def read_ome(path, level=0):
    """
    Input:
        path: Path to the OME-TIFF file
        level: Resolution level to read (default: 0, highest resolution)
    Outputs:
        arr: numpy array of the image data at the specified resolution level
        xml: OME-XML metadata as a string
        axes: axes string describing the order of dimensions in arr
    """
    with tifffile.TiffFile(path) as tif:
        xml = tif.ome_metadata
        series = tif.series[0]
        axes = series.axes
        if level == 0:
            arr = series.asarray()
        else:
            arr = series.levels[level].asarray()
    return arr, xml, axes

def spatial_axes(axes):
    """
    Return (y_axis,x_axis) index positions from a tiffile axes string.
    Example:
        'CYX' -> (1,2)
    """
    try:
        return axes.index('Y'), axes.index('X')
    except ValueError:
        raise ValueError(f"axes string {axes} does not contain both 'Y' and 'X'")

def channel_axis(axes):
    """
    Return the index of the channel axis from a tifffile axes string, or None
    Example:
        'YXS' -> 2, 'CYX' -> 0, 'YX' -> None
    """
    for a in ("S", "C"):
        if a in axes:
            return axes.index(a)
    return None

def ome_tostring(root):
    """
    Returns the OME-XML as a string with the default namespace set to the OME namespace.
    """
    ET.register_namespace("", ome_ns)
    return ET.tostring(root, encoding="unicode")

def single_channel_metadata(xml, ome_type="uint8", name=None):
    """
    Returns modified OME-XML (after stain deconvolution) as a string.
    """
    if xml is None:
        return None
    root = ET.fromstring(xml)
    ns = {"ome": ome_ns}
    for pixels in root.findall(".//ome:Pixels", ns):
        pixels.set("SizeC", "1") # greyscale/single-channel after deconvolution
        pixels.set("Type", ome_type)
        # RGB samples are interleaved, while a single channel is not
        pixels.attrib.pop("Interleaved", None)
        # Only keep the first channel element
        channels = pixels.findall("ome:Channel", ns)
        for extra in channels[1:]:
            pixels.remove(extra)
        if not channels:
            # Channel has to come first in the OME schema's element order
            channel = ET.Element(f"{{{ome_ns}}}Channel", {"ID": "Channel:0:0"})
            pixels.insert(0, channel)
            channels = [channel]
        channels[0].set("SamplesPerPixel", "1")
        # -1 is the OME default, opaque white, suitable for single grayscale channel
        channels[0].set("Color", "-1")
        if name is not None:
            channels[0].set("Name", name)
        # TiffData maps IFDs to planes and describes the old, multi-channel layout
        for tiffdata in pixels.findall("ome:TiffData", ns):
            pixels.remove(tiffdata)
    return ome_tostring(root)

def set_interleaved(xml, interleaved=True):
    """
    Set Pixels/@Interleaved to match how the samples are actually stored.
    """
    if xml is None:
        return None
    root = ET.fromstring(xml)
    ns = {"ome": ome_ns}
    for pixels in root.findall(".//ome:Pixels", ns):
        pixels.set("Interleaved", "true" if interleaved else "false")
        # TiffData maps IFDs to planes and describes the old, non-interleaved layout
        for tiffdata in pixels.findall("ome:TiffData", ns):
            pixels.remove(tiffdata)
    return ome_tostring(root)

def swap_xy_metadata(xml):
    """
    Swap the X and Y metadata in the OME-XML for all Pixels elements.
    Use when the image array has been transposed or rotated.
    """
    if xml is None:
        return None
    root = ET.fromstring(xml)
    ns = {"ome": ome_ns}
    for pixels in root.findall(".//ome:Pixels", ns):
        for x_attr, y_attr in xy_pairs:
            x_val = pixels.get(x_attr)
            y_val = pixels.get(y_attr)
            for attr, val in (y_attr, x_val), (x_attr, y_val):
                if val is not None:
                    pixels.set(attr, val)
                else:
                    pixels.attrib.pop(attr, None)
    return ome_tostring(root)

def plane_count(arr, photometric):
    """
    Number of IFDs a write of `arr` produces: everything above the two spatial
    axes, or the three of them when the samples are interleaved.
    """
    spatial = 3 if photometric == "rgb" else 2
    planes = 1
    for n in arr.shape[:-spatial]:
        planes *= n
    return planes

def single_file_metadata(xml, planes):
    """
    Replace TiffData that names external files with one entry describing the
    single file being written.

    The DAPI image in a Xenium morphology_set is a multi-file OME, 
    so single_file_metadata rewrites its OME-XML to describe only a single file.
    """
    if xml is None:
        return None
    root = ET.fromstring(xml)
    ns = {"ome": ome_ns}
    external = [
        pixels for pixels in root.findall(f".//{{{ome_ns}}}Pixels")
        if any(td.find("ome:UUID[@FileName]", ns) is not None
               for td in pixels.findall(f"{{{ome_ns}}}TiffData"))
    ]
    if not external:
        return xml
    for pixels in external:
        for tiffdata in pixels.findall(f"{{{ome_ns}}}TiffData"):
            pixels.remove(tiffdata)
        # TiffData comes after Channel in the OME schema's element order
        channels = pixels.findall(f"{{{ome_ns}}}Channel")
        at = list(pixels).index(channels[-1]) + 1 if channels else 0
        pixels.insert(at, ET.Element(
            f"{{{ome_ns}}}TiffData",
            {"IFD": "0", "FirstZ": "0", "FirstC": "0", "FirstT": "0",
             "PlaneCount": str(planes)},
        ))
    # This UUID identifies a source file, not the one being written
    root.attrib.pop("UUID", None)
    return ome_tostring(root)

def ascii_xml(xml):
    """
    Return xml with non-ASCII escaped as numeric character references
    """
    if xml is None:
        return None
    return xml.encode("ascii", "xmlcharrefreplace").decode("ascii")

def write_utf8_description(path, xml, escaped):
    """
    Rewrite the ImageDescription of an already-written file as raw UTF-8.
    """
    if xml is None or escaped is None or xml == escaped:
        return
    with tifffile.TiffFile(path, mode="r+") as tif:
        tif.pages[0].tags["ImageDescription"].overwrite(xml.encode("utf-8"))

def default_axes(arr, photometric):
    """
    tifffile axes string for an array written with photometric.
    'YX' for 2D, 'YXS' for interleaved samples, otherwise 'CYX'-style.
    """
    if photometric == "rgb":
        return "YXS" if arr.ndim == 3 else "Q" * (arr.ndim - 3) + "YXS"
    if arr.ndim == 2:
        return "YX"
    return "C" + "Q" * (arr.ndim - 3) + "YX"

def ome_metadata_dict(xml, axes, channel_color=True):
    """
    The parts of xml worth keeping, as a tifffile metadata dict.

    Passing this to tifffile as `metadata=` makes it generate the OME-XML
    itself, rather than carrying a source description through. make_pyramid.py
    writes that way and its files read correctly in QuPath, where a Xenium
    description passed through does not -- the source carries Plate, Instrument
    and StructuredAnnotations that a plain pass-through preserves and that
    something in the reader chain objects to. Dimensions, pixel size and channel
    names all survive; the acquisition extras do not.
    """
    meta = {"axes": axes}
    if xml is None:
        return meta
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return meta
    pixels = root.find(f".//{{{ome_ns}}}Pixels")
    if pixels is None:
        return meta
    for attr in ("PhysicalSizeX", "PhysicalSizeY"):
        try:
            meta[attr] = float(pixels.get(attr))
        except (TypeError, ValueError):
            meta.pop(attr, None)
            continue
        unit = pixels.get(attr + "Unit")
        if unit is not None:
            meta[attr + "Unit"] = unit
    channels = pixels.findall(f"{{{ome_ns}}}Channel")
    channel = {}
    names = [c.get("Name") for c in channels]
    if names and all(n is not None for n in names):
        channel["Name"] = names
    colors = [c.get("Color") for c in channels]
    # Keep Color when the source set it: single_channel_metadata writes -1 so
    # that VALIS finds a colormap entry instead of dropping channel colours.
    if channel_color and colors and all(c is not None for c in colors):
        channel["Color"] = [int(c) for c in colors]
    if channel:
        meta["Channel"] = channel
    return meta

def write_ome_tif(path, arr, xml, compression="LZW", photometric=None, tile=(1024,1024),
                  bigtiff=True, rebuild=False, axes=None):
    """
    Write an OME-TIFF with the given array and OME-XML metadata.

    rebuild=True hands tifffile a metadata dict and lets it generate the OME-XML,
    instead of passing xml through as the description. See ome_metadata_dict.
    """
    if photometric is None:
        # RGB if last dim has 3 or 4 channels
        photometric = "rgb" if (arr.ndim == 3 and arr.shape[-1] in (3,4)) else "minisblack"

    if rebuild:
        metadata = ome_metadata_dict(xml, axes or default_axes(arr, photometric))
        resolution = resolution_from_xml(xml)
        with tifffile.TiffWriter(path, bigtiff=bigtiff, ome=True) as tif:
            tif.write(
                arr,
                photometric=photometric,
                compression=compression,
                tile=tile,
                metadata=metadata,
                resolution=resolution,
                resolutionunit="CENTIMETER" if resolution else None,
            )
        return

    xml = single_file_metadata(xml, plane_count(arr, photometric))
    resolution = resolution_from_xml(xml)
    escaped = ascii_xml(xml)

    with tifffile.TiffWriter(path, bigtiff=bigtiff, shaped=False) as tif:
        tif.write(
            arr,
            photometric=photometric,
            compression=compression,
            tile=tile,
            description=escaped,
            resolution=resolution,
            resolutionunit="CENTIMETER" if resolution else None,
        )
    write_utf8_description(path, xml, escaped)

def pyramid_levels(arr, y_axis_idx, x_axis_idx, min_dim=256, max_levels=10):
    """
    Returns number of halvings needed to bring the smaller spatial 
    side down to min_dim.
    """
    h, w = arr.shape[y_axis_idx], arr.shape[x_axis_idx]
    n = 0
    while min(h, w) > min_dim and n < max_levels:
        h //= 2
        w //= 2
        n += 1
    return n

"""
def write_pyramidal_ome_tif(path, arr, xml, y_axis_idx, x_axis_idx, subresolutions=None, compression="zlib", photometric=None, tile=(1024, 1024), bigtiff=True, maxworkers=4):
    '''
    Write a pyramidal, tiled OME-TIFF (required by Xenium Explorer)
    Xenium Explorer accepts ZLIB or lossless JPEG 2000.

    '''
    if photometric is None:
        photometric = "rgb" if (arr.ndim == 3 and arr.shape[-1] in (3, 4)) else "minisblack"
    if subresolutions is None:
        subresolutions = pyramid_levels(arr, y_axis_idx, x_axis_idx)

    xml = single_file_metadata(xml, plane_count(arr, photometric))
    resolution = resolution_from_xml(xml)

    escaped = ascii_xml(xml)

    options = dict(photometric=photometric, compression=compression, tile=tile, maxworkers=maxworkers)
    if resolution:
        options["resolutionunit"] = "CENTIMETER"

    # Strided subsampling
    downsample = [slice(None)] * arr.ndim
    downsample[y_axis_idx] = slice(None, None, 2)
    downsample[x_axis_idx] = slice(None, None, 2)
    downsample = tuple(downsample)

    with tifffile.TiffWriter(path, bigtiff=bigtiff, shaped=False) as tif:
        tif.write(arr, subifds=subresolutions, description=escaped,
                  resolution=resolution, **options)
        scale = 1.0
        for _ in range(subresolutions):
            arr = arr[downsample]
            # Each level is half-scale, so its pixels cover twice the distance
            scale /= 2
            tif.write(arr, subfiletype=1,
                      resolution=tuple(r * scale for r in resolution) if resolution else None,
                      **options)
    write_utf8_description(path, xml, escaped)
"""

