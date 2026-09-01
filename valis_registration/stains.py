import numpy as np
from skimage.color import rgb2hed

def hematoxylin_deconvolution(arr):
    """
    Input:
        arr: 3D numpy array representing an RGB image (Y, X, 3)
    Output:
        2D numpy array representing the hematoxylin channel of the HED color space
    """
    print("Color deconvolution (RGB -> HED)")
    # Separate the stains from the H&E
    # rgb2hed output is 3D array with shape (Y, X, 3)
    hed = rgb2hed(arr) 
    hematoxylin = hed[:, :, 0]

    # Rescale hematoxylin to 0-255 uint8
    # Use percentile clipping to avoid extreme values
    lo, hi = np.percentile(hematoxylin, (1, 99))
    hematoxylin = np.clip(hematoxylin, lo, hi)
    hematoxylin = ((hematoxylin - lo) / (hi - lo) * 255).astype(np.uint8)

    return hematoxylin