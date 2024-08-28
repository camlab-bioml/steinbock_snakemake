import anndata as ad
import phenograph
import argparse
from pathlib import Path
# Parser object
parser = argparse.ArgumentParser(prog='Phenograph clustering',
                                 description='Perform phenograph clustering on an h5ad intensity object.')

# I/O Arguments
parser.add_argument("-i", "--input", dest="input", action="store", default="", help="Input h5ad file.")
parser.add_argument("-o", "--output", dest="output", action="store", default="", help="Output h5ad file.")
parser.add_argument("-k", "--nearest-neighbors", dest="k", action="store", type=int, default=30,
                    help="K value nearest neighbors to use in phenograph. Default is 30")
parser.add_argument("-m", "--min-cluster-size", dest="min_clust", action="store", type=int, default=10,
                    help="Minimum cluster size to use in phenograph. Default is 10")
parser.add_argument("-n", "--normalize", dest="normalize", action="store_true",
                    help="Column normalize the expression intensities prior to clustering.")

args = parser.parse_args()

export = ad.read_h5ad(args.input)

# column normalize?
expr = ((export.X - export.X.min()) / (export.X.max() - export.X.min())) if args.normalize else export.X

communities, graph, Q = phenograph.cluster(expr, k=args.k, min_cluster_size=args.min_clust)

export.obs['phenograph'] = communities

export.write_h5ad(Path(args.output))