import numpy as np
import skimage
import glob
import os
import skimage.io

input_folder = snakemake.input['raw_folder']
output_folder = snakemake.output["arcsinh_folder"]

# Create output director
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Get tiff filenames
filenames = glob.glob(input_folder + "/*.tiff")

# Read in tiffs and perform arcsinh normalization
for i in filenames:
    raw = skimage.io.imread(i)
    output_file = i.replace("raw", "arcsinh")
    output_file = output_file.replace(".tiff", "_arcsinh.tiff")
    normalized = np.arcsinh(raw)
    skimage.io.imsave(arr = normalized, fname = output_file)
    
exit()