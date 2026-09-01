"""
Use to make an image compatible with Xenium Explorer.
Converts a TIFF image into a pyramidal OME-TIFF with JPEG 2000 compression and 1024x1024 tiles.

Adapted from the code from:
https://www.10xgenomics.com/support/software/xenium-explorer/latest/tutorials/xe-image-file-conversion

Submit with make_pyramid.sh
"""

import tifffile as tf
import sys
import numpy as np
import cv2
def img_resize(img,scale_factor):
    width = int(np.floor(img.shape[1] * scale_factor))
    height = int(np.floor(img.shape[0] * scale_factor))
    return cv2.resize(img, (width, height), interpolation = cv2.INTER_AREA)
def write_ome_tif(filename, image, channel_names, photometric_interp, metadata, subresolutions):
    subresolutions = subresolutions
    fn = filename + "_xenium_compat.ome.tif"
    with tf.TiffWriter(fn,  bigtiff=True, ) as tif:
        px_size_x = metadata['PhysicalSizeX']
        px_size_y = metadata['PhysicalSizeY']
        options = dict(
            photometric=photometric_interp,
            tile=(1024, 1024),
            maxworkers=4,
            compression='jpeg2000',
            compressionargs={'level':85},
            resolutionunit='CENTIMETER',
        )
        print("Writing pyramid level 0")
        tif.write(
            image,
            subifds=subresolutions,
            resolution=(1e4 / px_size_x, 1e4 / px_size_y),
            metadata=metadata,
            **options
        )
        scale = 1
        for i in range(subresolutions):
            scale /= 2
            if photometric_interp == 'minisblack':
                if image.shape[0] < image.shape[-1]:
                    image = np.moveaxis(image,0,-1)
                    image = img_resize(image,0.5)
                    image = np.moveaxis(image,-1,0)
            else:
                image = img_resize(image,0.5)
            print("Writing pyramid level {}".format(i+1))
            tif.write(
                image,
                subfiletype=1,
                resolution=(1e4 * scale / px_size_x, 1e4 * scale / px_size_y),
                **options
            )
with tf.TiffFile(sys.argv[1]) as tif:
    image = tif.asarray()
    channel_names=None
    if tif.pages[0].samplesperpixel == 3:
        photometric_interp='rgb'
    elif tif.pages[0].samplesperpixel == 1:
        photometric_interp='minisblack'
    if tif.pages[0].get_resolution():
        res = tif.pages[0].get_resolution()
        if tif.pages[0].resolutionunit == 2:
            unit = "in"
        elif tif.pages[0].resolutionunit == 3:
            unit = "cm"
        elif tif.pages[0].resolutionunit == 4:
            unit = "mm"
        elif tif.pages[0].resolutionunit == 5:
            unit = "um"
    if tif.is_ome:
        meta = tf.xml2dict(tif.ome_metadata)
        res = (meta['OME']['Image']['Pixels']['PhysicalSizeX'],meta['OME']['Image']['Pixels']['PhysicalSizeY'])
        unit = meta['OME']['Image']['Pixels']['PhysicalSizeXUnit']
        try:
            channel_names=[]
            for i, element in enumerate(meta['OME']['Image']['Pixels']['Channel']):
                channel_names.append(meta['OME']['Image']['Pixels']['Channel'][i]['Name'])
        except KeyError:
            channel_names=None
    metadata={
        'PhysicalSizeX': res[0],
        'PhysicalSizeXUnit': unit,
        'PhysicalSizeY': res[1],
        'PhysicalSizeYUnit': unit,
        'Channel': {'Name': channel_names}
        }
filename = sys.argv[1]
if filename.endswith('.ome.tif'):
    filename = filename[:-len('.ome.tif')]
else:
    filename = filename.rsplit('.',1)[0]
write_ome_tif(filename, image, channel_names, photometric_interp, metadata, subresolutions = 7)