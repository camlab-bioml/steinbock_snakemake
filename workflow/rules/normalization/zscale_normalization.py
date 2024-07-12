import numpy as np
import skimage
import glob
import os
import cv2
import skimage.io

input_file = snakemake.input['raw']
output_file = snakemake.output["normalized"]

f = snakemake.params['alpha']
b = snakemake.params['beta']

# Read in tiffs and perform zscale normalization
raw = skimage.io.imread(input_file)
normalized = raw.copy()
for j in range(raw.shape[0]):
    normalized[j,:,:] = cv2.convertScaleAbs(raw[j,:,:], alpha = f, beta = b)
    
# Write tiffs to disk
skimage.io.imsave(arr = normalized, fname = output_file)

exit()