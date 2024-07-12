import skimage
import cv2
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cellSAM

img = "data/data_analysis_training/img/zscale/Melanoma_1um_001.tiff"
panel = "data/data_analysis_training/panel.csv"

nuclear = ["Ir191"]
membrane = ["Er170"]

img = skimage.io.imread(img)
panel = pd.read_csv(panel)

# Load in nuclear channel
nuclear_index = panel[panel['channel'] == 'Ir193'].index
nuclear_tiff = img[nuclear_index,:,:]

# Load in membrane channels
membrane_index = panel[panel['channel'].isin(membrane)].index
membrane_tiff = img[membrane_index,:,:]
membrane_average = np.mean(membrane_tiff, axis=0)

# Merge together to form rgb image, nuclear = blue, membrane = green, 
merged_image = np.zeros((nuclear_tiff.shape[1], nuclear_tiff.shape[2], 3), dtype=np.uint64)
merged_image[:, :, 1] = membrane_average
merged_image[:, :, 2] = nuclear_tiff

# Tile
tile_size = (25, 25)

mask, embedding, bounbding_boxes= cellSAM.segment_cellular_image(merged_image, device='cuda', normalize=False)

plt.imshow(merged_image)
plt.show()
plt.close()