import os
import matplotlib.pyplot as plt
import skimage.io
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

raw_img = snakemake.input['raw']
arcsinh_img = snakemake.input['arcsinh']
zscale_img = snakemake.input['zscale']
panel = snakemake.input['panel']
projects = snakemake.wildcards['projects']
roi = snakemake.wildcards['ROI']

path = "data/" + projects + "/img/plots/" + roi
os.makedirs(path, exist_ok=True)

# Read in tiffs and channel metadata
raw_img = skimage.io.imread(raw_img)
arcsinh_img = skimage.io.imread(arcsinh_img)
zscale_img = skimage.io.imread(zscale_img)

panel = pd.read_csv(panel)
panel = panel['name'].to_list()

for i in range(0, raw_img.shape[0]):
    data = [raw_img[i,:,:], arcsinh_img[i,:,:], zscale_img[i,:,:]]
    title = ["raw - " + panel[i], "arcsinh - " + panel[i], "zscale - " + panel[i]]
    filename = path + "/" + roi + "_" + panel[i] + ".png"
    fig, axs = plt.subplots(1, 3, figsize=(30, 8), layout="constrained")
    # Iterate over each subplot, data, colormap, and title
    for j, ax in enumerate(axs):
        cax = ax.imshow(data[j], cmap="gray")
        divider = make_axes_locatable(ax)
        cax_divider = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(cax, cax = cax_divider)
        ax.set_title(title[j])
    plt.savefig(filename, dpi = 300)
    plt.close()

exit()