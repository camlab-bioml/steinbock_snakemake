import os
import pandas as pd
from PIL import Image
import numpy as np
import anndata as ad
import unittest
import pytest
import glob
import tifffile
import json


class SteinbockSnakemakeIntegrationTests(unittest.TestCase):

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, get_steinbock_out_dir):
        self.get_steinbock_out_dir = get_steinbock_out_dir

    @pytest.mark.usefixtures("get_steinbock_out_dir")
    def test_check_panel(self):
        panel = pd.read_csv(os.path.join(self.get_steinbock_out_dir, 'panel.csv'))
        assert panel.shape == (12, 6)

    @pytest.mark.usefixtures("get_steinbock_out_dir")
    def test_check_images(self):
        image_dir = os.path.join(self.get_steinbock_out_dir, 'img', 'raw')
        raw = tifffile.imread(glob.glob(f'{image_dir}/*.tiff')[0])
        assert raw.shape == (12, 500, 500)


    @pytest.mark.usefixtures("get_steinbock_out_dir")
    def test_check_masks(self):
        mask = np.array(Image.open(os.path.join(self.get_steinbock_out_dir,
            'deepcell', 'cell', 'test_018.tiff'))).astype(np.uint32)
        assert np.max(mask) == 581

    @pytest.mark.usefixtures("get_steinbock_out_dir")
    def test_check_intensities(self):
        intensities = pd.read_csv(os.path.join(self.get_steinbock_out_dir, 'quantification',
                                           'intensities', 'test_018.csv'))

        assert intensities.shape == (581, 13)

    @pytest.mark.usefixtures("get_steinbock_out_dir")
    def test_check_export(self):
        export_anndata = ad.read_h5ad(os.path.join(self.get_steinbock_out_dir,
                                               'export', 'test_mcd.h5ad'))

        assert export_anndata.X.shape == (581, 12)

        assert 'phenograph' in export_anndata.obs

        assert 'UMAP' in export_anndata.obsm
        assert export_anndata.obsm['UMAP'].shape == (581, 2)

        umap_coordinates = pd.read_csv(os.path.join(self.get_steinbock_out_dir, 'export',
                                           'umap_coordinates.csv'))
        assert umap_coordinates.shape == (581, 2)

    @pytest.mark.usefixtures("get_steinbock_out_dir")
    def test_roi_scaling(self):
        with open(os.path.join(self.get_steinbock_out_dir,
                               'export', 'scaling.json')) as json_file:
            scaling = json.load(json_file)
            assert isinstance(scaling, dict)
            assert 'channels' in list(scaling.keys())
            assert len(scaling['channels']) == 12
