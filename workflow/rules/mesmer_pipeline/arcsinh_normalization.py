import numpy as np
import skimage
import glob
import os
import skimage.io

input_folder = snakemake.input[0] + "/" 
output_folder = snakemake.output["normalized"] + "/" 

# Get Tiff files from steinbock output
files = glob.glob(input_folder + '*.tiff')

# Generate output filenames
filenames = []
for i in range(len(files)):
    new_name = output_folder + "/" + os.path.basename(files[i])
    filenames.append(new_name)

# Read in tiffs and perform arcsinh normalization
tiffs = []
for i in files:
    array = skimage.io.imread(i)
    for j in range(array.shape[0]):
            array[j,:,:] = np.arcsinh(array[j,:,:])
    tiffs.append(array)
    
# make output folder
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
# Write tiffs to disk
for i in range(len(tiffs)):
    skimage.io.imsave(arr = tiffs[i], fname = filenames[i])

exit()