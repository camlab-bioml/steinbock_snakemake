# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.0.6]

### Added

- `mask_id` segmentation ID column now part of `export.h5ad` .obs
- RGB of nuclear markers + segmentation overlay in `deepcell`, `overlay` #28

### Changed

- Conditional parse of txt files or mcd (mcd by default) using config `process_txt` #26
- scaling JSON for rakaia now supports parsing .txt files #33

### Fixed

- raw tiff ROIs that trigger a `panopticnet` `ValueError` on cell segmentation 
will be moved to `raw_not_quantified` #36

### Removed

- Duplicate script files and directories

## [v0.0.5]

### Added
- Updated environment.yaml file for quick conda install.
- Likely frozen version for main.
- Snakemake profiles for 1 line pipeline execution with configured snakemake options for newer versions of snakemake >= 8. 

### Changed
- Verbosity improvements for snakemake rule pointers

### Removed
- DAG Removed (Not much purpose)
- archsinh, zscale and mixmax previews removed for lightweight pipeline. 
- Removed srcdir() functions in quantification rules as they are deprecated.

## [v0.0.4]

### Added

- JSON output for auto-scaling based on all ROIs for rakaia visualization
- multi UMAP coordinate output for minimum distance + plotting w/ phenograph
- Input channels to ignore for UMAP and phenograph clustering #22

## [v0.0.3]

### Added

- clustering using [phenograph](https://github.com/dpeerlab/PhenoGraph) (added to exported h5ad)
- basic UMAP projection (coordinates output in CSV)

## [v0.0.2]

### Added

- Simplified Snakemake rule structure
- Zarr store as well as h5ad outputs

### Changed

 - From container based workflow to conda and pip based workflow. To deal with docker container issues.

### Removed

 - Requirement of containers


## [v0.0.1]

### Added

- DeepCell Segmentation using mesmer
- Pipeline created

### Changed

 - NULL

### Removed

 - NULL