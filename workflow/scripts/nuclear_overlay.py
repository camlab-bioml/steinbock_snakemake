"""
Generate a png overlay of every ROI's nuclear channels with the mask outlines
"""

import glob
import warnings
import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from PIL import Image
from tifffile import TiffFile, imwrite
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay


def cli_parser():
    parser = argparse.ArgumentParser(add_help=False,
            description="Generate a png overlay of every ROI's nuclear channels with the mask outlines",
            usage='Example:\n nuclear_overlay.py -i raw_tiff_directory_input -o deepcell/overlay')
    parser.add_argument('-i', "--input", action="store",
                        help="Path to raw img input (one tiff per ROI)",
                        dest="input", type=str, required=True)
    parser.add_argument('-m', "--mask", action="store",
                        help="Path to segmentation cell mask files (one tiff per ROI)",
                        dest="mask", type=str, required=True)
    parser.add_argument('-o', "--outdir", action="store",
                        help="Set the output tiff file. Default is geojson.tiff written to the current directory",
                        dest="outdir", default="deepcell/overlay", type=str)
    parser.add_argument("-p", "--panel", dest="panel", action="store", default="data/data_analysis_training/panel.csv",
                        help="Name of the deepcell panel (should be panel_deepcell.csv)")
    parser.add_argument('-t', "--type", action="store",
                        help="Set the type of image output (either png or tiff)",
                        dest="type", type=str, required=True, default="png")

    return parser

def main(sysargs=sys.argv[1:]):
    warnings.filterwarnings("ignore")
    parser = cli_parser()
    args = parser.parse_args(sysargs)

    panel_dp = pd.read_csv(args.panel)
    nuclear_indices = panel_dp[panel_dp['deepcell'] == 1].index.to_list()
    raw_tiffs = glob.glob(f"{args.input}/*.tiff")
    cols_channel = ['green'] * len(nuclear_indices)
    for roi in raw_tiffs:
        # Make sure to get the basename to get the matching mask
        roi_base = Path(roi).stem
        with TiffFile(roi) as tiff_open:
            img_array = tiff_open.asarray()
            C, H, W = img_array.shape

            only_nuc_channels = img_array[nuclear_indices, :, :]
            only_nuc_channels = np.stack([arr for arr in only_nuc_channels], axis=2)

            rgb_nuc = create_rgb_image(only_nuc_channels.reshape(1, H, W, len(nuclear_indices)),
                                       channel_colors=cols_channel)

        mask = f"{args.mask}/{roi_base}.tiff"
        with TiffFile(mask) as mask_open:
            mask = mask_open.asarray().astype(np.uint32)
            H_mask, W_mask = mask.shape

        overlay_data = make_outline_overlay(rgb_data=rgb_nuc,
                    predictions=mask.reshape(1, H_mask, W_mask, 1))

        out_img = Image.fromarray((overlay_data[0] * 255).clip(0, 255).astype(np.uint8)).convert('RGB')
        if args.type == 'tiff':
            imwrite(f"{args.outdir}/{roi_base}.tiff",
                    np.array(out_img).astype(np.uint8), photometric='RGB')
        else:
            out_img.save(f"{args.outdir}/{roi_base}.png", format="PNG")

if __name__ == "__main__":
    main()