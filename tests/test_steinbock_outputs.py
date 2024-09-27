import os
import pandas as pd
from PIL import Image
import numpy as np
import anndata as ad

def test_steinbock_outputs(get_steinbock_out_dir):
    panel = pd.read_csv(os.path.join(get_steinbock_out_dir, 'panel.csv'))
    assert panel.shape == (12, 6)

    mask = np.array(Image.open(os.path.join(get_steinbock_out_dir,
            'deepcell', 'cell', 'test_018.tiff'))).astype(np.uint32)
    assert np.max(mask) == 581

    intensities = pd.read_csv(os.path.join(get_steinbock_out_dir, 'quantification',
                                           'intensities', 'test_018.csv'))

    assert intensities.shape == (581, 13)

    export_anndata = ad.read_h5ad(os.path.join(get_steinbock_out_dir,
                                               'export', 'test_mcd.h5ad'))
    assert export_anndata.X.shape == (581, 12)

    assert 'phenograph' in export_anndata.obs

    umap_coords = pd.read_csv(os.path.join(get_steinbock_out_dir, 'export',
                                           'umap_coordinates.csv'))
    assert umap_coords.shape == (581, 2)
    