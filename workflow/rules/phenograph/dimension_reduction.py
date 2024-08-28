import anndata as ad
import argparse
import umap
from sklearn.preprocessing import StandardScaler
import pandas as pd
# Parser object
parser = argparse.ArgumentParser(prog='UMAP',
                                 description='Generate UMAP coordinates for objects from an intensity h5ad.')

# I/O Arguments
parser.add_argument("-i", "--input", dest="input", action="store", default="", help="Input h5ad file.")
parser.add_argument("-o", "--output", dest="output", action="store",
                    default="umap_coordinates.csv", help="Output CSV file containing UMAP1 and UMAP1 coordinates.")
parser.add_argument("-n", "--normalize", dest="normalize", action="store_true",
                    help="Column normalize the expression intensities prior to clustering.")

args = parser.parse_args()

export = ad.read_h5ad(args.input)

# column normalize?
expr = ((export.X - export.X.min()) / (export.X.max() - export.X.min())) if args.normalize else export.X

try:
    umap_obj = umap.UMAP()
except (AttributeError, ImportError):
    import umap
    umap_obj = umap.UMAP()
scaled = StandardScaler().fit_transform(expr)
embedding = pd.DataFrame(umap_obj.fit_transform(scaled), columns=["UMAP1", "UMAP2"])

embedding.to_csv(args.output, index=False)



