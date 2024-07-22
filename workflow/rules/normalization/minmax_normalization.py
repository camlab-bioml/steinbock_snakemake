import numpy as np
import glob
import os
import skimage

input_folder = snakemake.input['raw_folder']
output_folder = snakemake.output["minmax_folder"]

# Create output director
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Get tiff filenames
filenames = glob.glob(input_folder + "/*.tiff")

# Read in tiffs and perform minmax normalization
for i in filenames:
    raw = skimage.io.imread(i)
    output_file = i.replace("raw", "minmax")
    output_file = output_file.replace(".tiff", "_minmax.tiff")
    normalized = raw.copy()
    for j in range(raw.shape[0]):
        min = np.min(normalized[j,:,:])
        max = np.max(normalized[j,:,:])
        normalized[j,:,:] = (normalized[j,:,:] - min) / (max - min)
    skimage.io.imsave(arr = normalized, fname = output_file)
    
exit()