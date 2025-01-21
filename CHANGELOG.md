# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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