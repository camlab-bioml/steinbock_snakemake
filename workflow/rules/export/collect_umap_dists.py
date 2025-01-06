"""
Collect all the UMAP min distance coordinate lists into the export anndata
"""

import argparse


from pathlib import Path
import anndata as ad
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(prog='Collect UMAP distances',
        description='collect UMAP min distance coordinates into the `obsm` Anndata slot for the export h5ad.')

parser.add_argument("-i", "--input", dest="input", action="store", default="", help="Input h5ad file.")
parser.add_argument("-cd", "--coord-directory", dest="coord_dir", action="store",
                    help="Import directory for UMAP min distance coordinate lists.")
parser.add_argument("-o", "--output", dest="output", action="store",
                    help="Output anndata file with every UMAP min distance computation saved in `obsm`")

args = parser.parse_args()


coord_list = [str(i) for i in Path(args.coord_dir).rglob("*coordinates.csv")]

export = ad.read_h5ad(args.input)

for min_dist_coordinates in coord_list:
    idx1 = min_dist_coordinates.index('umap_')
    idx2 = min_dist_coordinates.index('_coordinates.csv')
    dist_val = min_dist_coordinates[idx1 + len('umap_') + 0: idx2]
    embedding = pd.read_csv(min_dist_coordinates)
    export.obsm[f'UMAP_{dist_val}'] = np.array(embedding).astype(np.float32)

export.write_h5ad(args.output)

