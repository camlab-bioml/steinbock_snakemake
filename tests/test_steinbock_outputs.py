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
from pathlib import Path


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
        assert 500 <= np.max(mask) <= 600

    @pytest.mark.usefixtures("get_steinbock_out_dir")
    def test_check_intensities(self):
        intensities = pd.read_csv(os.path.join(self.get_steinbock_out_dir, 'quantification',
                                           'intensities', 'test_018.csv'))

        assert intensities.shape[1] == 13
        assert 500 <= intensities.shape[0] <= 600

    @pytest.mark.usefixtures("get_steinbock_out_dir")
    def test_check_export(self):
        export_anndata = ad.read_h5ad(os.path.join(self.get_steinbock_out_dir,
                                               'export', 'test_mcd.h5ad'))

        assert 500 <= export_anndata.X.shape[0] <= 600
        assert export_anndata.X.shape[1] == 12

        assert 'phenograph' in export_anndata.obs

        for min_dist in [0, 0.1, 0.25, 0.5, 1]:
            assert f'UMAP_min_dist_{min_dist}' in export_anndata.obsm
            assert export_anndata.obsm[f'UMAP_min_dist_{min_dist}'].shape[1] == 2
            assert 500 <= export_anndata.obsm[f'UMAP_min_dist_{min_dist}'].shape[0] <= 600

        umap_coord_list = sorted([str(i) for i in Path(
            os.path.join(self.get_steinbock_out_dir, 'export')).rglob('*coordinates.csv')])

        assert len(umap_coord_list) == 5

        for umap_dist in umap_coord_list:
            umap_coordinates = pd.read_csv(umap_dist)
            assert 500 <= umap_coordinates.shape[0] <= 600
            assert umap_coordinates.shape[1] == 2

        umap_plot_list = sorted([str(i) for i in Path(
            os.path.join(self.get_steinbock_out_dir, 'export')).rglob('*.png')])

        for min_dist in [0, 0.1, 0.25, 0.5, 1]:
            assert (any([str(min_dist) in plot_file for plot_file in umap_plot_list]))

        assert len(umap_plot_list) == 5

    @pytest.mark.usefixtures("get_steinbock_out_dir")
    def test_roi_scaling(self):
        with open(os.path.join(self.get_steinbock_out_dir,
                               'export', 'scaling.json')) as json_file:
            scaling = json.load(json_file)
            assert isinstance(scaling, dict)
            assert 'channels' in list(scaling.keys())
            assert len(scaling['channels']) == 12
