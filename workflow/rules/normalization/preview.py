import os
import glob
import matplotlib.pyplot as plt
import skimage.io
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

raw_folder = snakemake.input['raw']
arcsinh_folder = snakemake.input['arcsinh']
zscale_folder = snakemake.input['zscale']
minmax_folder = snakemake.input['minmax']
panel = snakemake.input['panel']
images = snakemake.input['images']
projects = snakemake.wildcards['projects']
output_folder = "data/" + projects  + "/img/preview/"

# Read in tiffs and channel metadata
raw_files = glob.glob(raw_folder + f"/*.tiff")
arcsinh_files = glob.glob(arcsinh_folder + f"/*.tiff")
zscale_files = glob.glob(zscale_folder + f"/*.tiff")
minmax_files = glob.glob(minmax_folder + f"/*.tiff")

print(panel)
panel = pd.read_csv(panel)
panel = panel['name'].to_list()

images = pd.read_csv(images).image.to_list()
images = [s.replace(".tiff","") for s in images]

suffix = output_folder
output_folders = [suffix + s for s in images]

for i in output_folders:
    if not os.path.exists(i):
        os.makedirs(i)
            
for i in range(0,len(raw_files)):
    raw_img = skimage.io.imread(raw_files[i])
    minmax_img = skimage.io.imread(minmax_files[i])
    arcsinh_img = skimage.io.imread(arcsinh_files[i])
    zscale_img = skimage.io.imread(zscale_files[i])
    for j in range(0, raw_img.shape[0]):
        sub_plots = [raw_img[j,:,:], minmax_img[j,:,:], arcsinh_img[j,:,:], zscale_img[j,:,:]]
        titles = ['Raw ' + panel[j], 'minmax ' + panel[j], 'arcsinh ' + panel[j], 'zscale ' + panel[j]]
        filename = output_folders[i] + "/" + panel[j] + ".png"
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))
        for k, ax in enumerate(axes.flat):
            img = ax.imshow(sub_plots[k], cmap='gray')
            ax.set_title(titles[k])
            ax.set_xticks([])  # Remove x-axis ticks
            ax.set_yticks([])  # Remove y-axis ticks
    # Create a colorbar that spans the height of the plot
            cbar = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=10)
        plt.tight_layout()
        plt.savefig(filename, bbox_inches = 'tight')
        plt.close()
        



exit()