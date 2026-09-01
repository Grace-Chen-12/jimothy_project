import numpy as np

rigid = {
    "identity": lambda a, y, x: a,
    "rot90": lambda a, y, x: np.rot90(a, k=1, axes=(y, x)),
    "rot180": lambda a, y, x: np.rot90(a, k=2, axes=(y, x)),
    "rot270": lambda a, y, x: np.rot90(a, k=3, axes=(y, x)),
    "flipud": lambda a, y, x: np.flip(a, axis=y),
    "fliplr": lambda a, y, x: np.flip(a, axis=x),
    "transpose": lambda a, y, x: np.rot90(np.flip(a, axis=x), k=1, axes=(y, x)),
    "transverse": lambda a, y, x: np.rot90(np.flip(a, axis=y), k=1, axes=(y, x)),
}
swaps_axes_transforms = {"rot90", "rot270", "transpose", "transverse"}

def apply_transform(arr, name, y_axis_idx, x_axis_idx):
    if name not in rigid:
        raise ValueError(f"Unknown transform {name}. Must be one of: {list(rigid.keys())}")
    return rigid[name](arr, y_axis_idx, x_axis_idx)

    
def swaps_pixels(name):
    """True if transform exchanges X and Y axes"""
    return name in swaps_axes_transforms


