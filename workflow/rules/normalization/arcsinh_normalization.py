import numpy as np
import skimage
import glob
import os
import skimage.io

input_file = snakemake.input['raw']
output_file = snakemake.output["normalized"]

# Read in tiffs and perform arcsinh normalization
raw = skimage.io.imread(input_file)
normalized = raw.copy()
for j in range(raw.shape[0]):
    normalized[j,:,:] = np.arcsinh(raw[j,:,:])
    
# Write tiffs to disk
skimage.io.imsave(arr = normalized, fname = output_file)

exit()