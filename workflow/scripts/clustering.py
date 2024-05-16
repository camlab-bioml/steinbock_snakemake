import spatialdata as sd
import spatialdata_io as sdio
import spatialdata_plot as sdplt
import os

zarr_store = "data/test_mcd/export/"

zarr = sdio.steinbock(zarr_store, labels_kind="deepcell", )