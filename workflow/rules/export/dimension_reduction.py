import anndata as ad
import argparse
import umap
from sklearn.preprocessing import StandardScaler
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

def set_plot_color_mapping(export_ad: ad.AnnData,
                           cluster_cat: str='phenograph'):

    if cluster_cat in list(export.obs.columns):

        export_ad.obs[cluster_cat] = export_ad.obs[cluster_cat].apply(str)

        # Get unique categories from the 'category' column
        unique_cats = export_ad.obs[cluster_cat].unique()

        # Assign colors using the tab20 colormap
        color_map = matplotlib.colormaps['tab20']
        cat_colors = {category: color_map(i / len(unique_cats)) for
                      i, category in enumerate(unique_cats)}

        # Map colors to each row based on the category
        mapped_cols = export_ad.obs[cluster_cat].map(cat_colors)

    else:
        mapped_cols = None
        cat_colors = {}

    label = None if mapped_cols is None else 'phenograph'
    return mapped_cols, cat_colors, label

parser = argparse.ArgumentParser(prog='UMAP',
                                 description='Generate UMAP coordinates for objects from an intensity h5ad.')

parser.add_argument("-i", "--input", dest="input", action="store", default="", help="Input h5ad file.")
parser.add_argument("-oc", "--output-coords", dest="output_coords", action="store",
                    default="umap_coordinates.csv", help="Output CSV file containing UMAP1 and UMAP1 coordinates.")
parser.add_argument("-n", "--normalize", dest="normalize", action="store_true",
                    help="Column normalize the expression intensities prior to clustering.")
parser.add_argument("-p", "--plot", dest="plot", action="store_true",
                    help="Plot the UMAP coordinates with export clustering overlaid.")
parser.add_argument("-op", "--output-plot", dest="output_plot", action="store",
                    help="Output png file for the scatter plot with export clustering.")
parser.add_argument("-md", "--min-dist", dest="min_dist", action="store", default=0.1,
                    type=float, help="Min dist value for UMAP")


args = parser.parse_args()

export = ad.read_h5ad(args.input)

# column normalize?
expr = ((export.X - export.X.min()) / (export.X.max() - export.X.min())) if args.normalize else export.X

try:
    umap_obj = umap.UMAP(min_dist=args.min_dist)
except (AttributeError, ImportError):
    import umap
    umap_obj = umap.UMAP(min_dist=args.min_dist)
scaled = StandardScaler().fit_transform(expr)
embedding = pd.DataFrame(umap_obj.fit_transform(scaled), columns=["UMAP1", "UMAP2"])

embedding.to_csv(args.output_coords, index=False)

if args.plot:

    mapping, category_colors, lab = set_plot_color_mapping(export)

    plt.figure(figsize=(16, 10))
    scatter = plt.scatter(
        embedding['UMAP1'],
        embedding['UMAP2'],
        c=mapping,
        s=12)

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title(f"UMAP Scatter Plot, {args.min_dist}")

    handles = [plt.Line2D([0], [0], marker='o', color=color, linestyle='', markersize=10)
               for category, color in category_colors.items()]
    plt.legend(handles, category_colors.keys(), title=lab, loc="best")

    # Save the plot
    plt.tight_layout()
    plt.savefig(args.output_plot)