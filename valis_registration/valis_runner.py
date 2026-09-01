import logging
import time
import os
from pathlib import Path

import io_utils
from valis import registration


class _TiffTagOrderFilter(logging.Filter):
    """
    Filter that removes libvips warnings about TIFF tag order.
    """
    def filter(self, record):
        return "tags are not sorted in ascending order" not in record.getMessage()

logging.getLogger("pyvips").addFilter(_TiffTagOrderFilter())


def stage_images(moving_arr, moving_xml, moving_name, reference_name, reference_path, staged_dir,
                 compression="LZW"):
    """Stages the input images for VALIS registration"""
    moving_out = os.path.join(staged_dir, moving_name)
    reference_out = os.path.join(staged_dir, reference_name)
    io_utils.write_ome_tif(moving_out, moving_arr, moving_xml, compression=compression)

    reference_arr, reference_xml, _ = io_utils.read_ome(reference_path)
    io_utils.write_ome_tif(reference_out, reference_arr, reference_xml, compression=compression)

def run_registration(moving, fixed, output_dir, kill_jvm=True):
    """
    Run VALIS registration
    """
    fixed_image_dir = os.path.dirname(fixed)
    registered_slides_dir = os.path.join(output_dir, "registered_slides")
    os.makedirs(registered_slides_dir, exist_ok=True)


    start = time.time()
    registrar = registration.Valis(
        # note src_dir functionally doesn't matter since the actual images are
        # specified in img_list, but it is required by the VALIS constructor.
        src_dir=fixed_image_dir, 
        dst_dir=output_dir,
        name="valis_results",
        img_list=[moving, fixed],
        reference_img_f=fixed,
        align_to_reference=True,
        crop="reference",
        # max_image_dim_px needs to be greater than max_processed_image_dim_px
        max_processed_image_dim_px=7999, 
        max_image_dim_px=8000, 
    )

    rigid_registrar, non_rigid_registrar, error_df = registrar.register()
    registrar.register_micro(max_non_rigid_registration_dim_px=8000)

    stop = time.time()
    elapsed = stop - start
    print(f"Registration took {elapsed/60:.2f} minutes")

    #  --- Guard against failed registration ---
    if error_df is None:
        print("Registration failed - error_df is None. Check logs above.")
        registration.kill_jvm()
        exit(1)
    print("\n--- Registration error summary ---")
    print(error_df.to_string(index=False))
    print()

    # --- Save registered images ---
    start = time.time()
    slide = registrar.get_slide(moving)
    slide.reader.max_image_dim_px = 2500
    registered_path = os.path.join(registered_slides_dir, f"{os.path.basename(moving).split('.')[0]}_registered.ome.tiff")
    try: 
        slide.warp_and_save_slide(registered_path, non_rigid=True, crop="reference",tile_wh=1024)
    except Exception as e:
        print(f"Failed to save {registered_path}: {e}")
        registration.kill_jvm()
        exit(1)
    stop = time.time()
    elapsed = stop - start
    print(f"Saving {registered_path} took {elapsed/60:.2f} minutes")

    # --- Shutdown the JVM ---
    if kill_jvm:
        registration.kill_jvm()

    return registered_path

def find_registrar_pickle(valis_results_dir):
    """
    Returns the path to the registrar pickle file.
    """
    data_dir = Path(valis_results_dir) / "data"
    # VALIS saves the registrar to <results_dir>/data/<name>_registrar.pickle.
    matches = sorted(data_dir.glob("*_registrar.pickle"))
    if not matches:
        raise FileNotFoundError(f"no *_registrar.pickle found in {data_dir}")
    if len(matches) > 1:
        raise ValueError(f"expected one registrar pickle in {data_dir}, found {[m.name for m in matches]}")
    return matches[0]

def register_pickle(was_registered, to_register, pickle_path, output_path):
    """
    Apply the transforms VALIS already learned for `was_registered` to a
    different file, `to_register`.
    """
    registrar = registration.load_registrar(str(pickle_path))
    slide = registrar.get_slide(str(was_registered))
    try:
        slide.warp_and_save_slide(
            str(output_path),
            src_f=str(to_register),
            non_rigid=True,
            crop="reference",
            tile_wh=1024,
        )
        print(f"Saved registered RGB H&E: {output_path}")
    except Exception as e:
        print(f"Failed to save {output_path}: {e}")
        registration.kill_jvm()
        raise
    registration.kill_jvm()
    return output_path

    