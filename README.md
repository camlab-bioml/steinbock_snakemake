# Steinbock-snakemake v0.0.1

## Pipeline overview
This steinbock pipeline uses mesmer to segment cells and nuclei, generate neighbors and finally outputs files required for downstream analysis. 

![DAG](dag.png)

## Set-up pipeline environment (via conda, pipenv)

We will first set up the pipeline environment via conda. This will install all dependencies including snakemake and singularity for use.

```
git clone https://github.com/ttj131/steinbock_snakemake && cd steinbock_snakemake # clone git repo, then move into the pipeline folder. 
conda env create -y --file=workflow/env/environment.yml
conda activate steinbock-snakemake
```

Alternatively we can use pipenv to just install snakemake given support for singularity in a cluster environment. (Not recommended)
```
pipenv install -r workflow/env/requirements.txt 
```

## Clone docker containers

### Cloning CPU-based container (Recommended)
After successful creation of the environment, we will then generate a singularity container for use with the snakemake pipeline. This container will require about ~4GB of storage, and will be located under workflow/envs

```
cd workflow/env
singularity pull steinbock-cpu.sif docker://ghcr.io/bodenmillergroup/steinbock:0.16.1
cd ../..
```

### Cloning GPU-based container (Not Recommended)
This will clone the GPU container for steinbock. Note, you will require `nvidia-container-toolkit` and require CUDA drivers (CUDA 11.8) to run the pipeline. Only use this for large datasets that require GPU support.
```
cd workflow/env
singularity pull steinbock-gpu.sif docker://ghcr.io/bodenmillergroup/steinbock:0.16.1-gpu
cd ../..
```

## Running Steinbock-snakemake

### Data input and output structures
Now that the containers and environment is set-up, we will run and test the pipeline on a small mcd file. This is the structure of the input folder:
```
── 📁data
    └── 📁test_mcd
        └── 📁mcd # Folder for you XTi mcd file or Hyperion mcd file. .txt files are also supported.
            └── test.mcd
```

### Creating a new project
Each project is given a directory under the data folder. To create a new project, make a new folder under `data/yourproject` and add a `data/yourproject/mcd` folder. Drag and drop your `yourproject.mcd` into the mcd folder. The minimum requirement for the project folder structure is as shown:
```
── 📁data
    └── 📁yourproject
        └── 📁mcd # Folder for you XTi *.mcd file or Hyperion *.mcd file. Hyperion *.txt files are also supported.
            └── test.mcd
```

### Preflight configuration
Before running the pipeline, we will first take a look at the configuration file `config/config.yaml`. This will store important information regarding the settings for each analysis module. Note that values for Rule Group 2 and Rule Group 3 are commented out, as we will not be using CellProfiler or Ilastik for this tutorial. The structure of the config file is as follows:

```
projects:
  - test_mcd # Informs the snakemake pipeline which projects to process

# Rule Group 1 : Options for preprocessing IMC Files
hpf: 5 # hot pixel filter

# Rule Group 4 : Deepcell
cytoplasm: "Gd158"
deepcell_app: "mesmer"
deepcell_model: "MultiplexSegmentation"
deepcell_modelpath: "/opt/keras/models"
deepcell_pxsize: 1

# Rule Group 5 :
aggr: "median"
neighbor_type: "borders"
dmax: 15
kmax: 5
```

Under the `projects` configuration, add your folder name and the pipeline will process the mcd files in that folder, if it had not already done so. We will use the `test_mcd` project in this demonstration. 
```
projects:
    - test_mcd
    # - myproject # <-- add the folder name to your project
```

### Channels for segmenting cells
To help mesmer segment cells, steinbock will require the user to input the mass channels for cytoplasm specific markers (e.g panCK, Actin etc). You should include those channels under `cytoplasm` in Rule Group 4.
```
cytoplasm: "Gd158" # <-- add mass channels that corresspond to the cytoplasm
deepcell_app: "mesmer"
deepcell_model: "MultiplexSegmentation"
deepcell_modelpath: "/opt/keras/models"
deepcell_pxsize: 1 # <-- pixel sizes for deepcell, default 1 um for IMC images. 
```

### Extracting features
After generating segmenting masks, steinbock will proceeed to extract single-cell features (intensities, neighbors, region properties etc.) The following options are avaialble under Rule Group 5
```
aggr: "median" # [sum|min|max|mean|median|std|var] are available options. We will be using the median intensity per cell. 
neighbor_type: "borders" # [centroids|borders|expansion] # Method for determining neighbors. By default determines thresholding on Euclidean distances between object borders.
dmax: 15 # Max distance
kmax: 5 # Max number of neighbors
```

### Project specific configuration
After editing `config.yml`, we would want to ensure we use the sames settings for each project. A simple solution is to copy the `config.yml` file over to the project directory before calling the snakemake, and directing it to the config file. 
```
cp config/config.yaml data/test_mcd/test_mcd.yaml
```

### Running the snakemake pipeline
After completing the preflight configuration, we can run the pipeline. The `snakemake` command is invoked in the `steinbock-snakemake` directory:
1. Specify cores required using the `-c` flag.
2. Specify the configuration file used using the `--configfile` flag.
3. Specify bind mount for singularity container using the `--singularity-args` flag. This is required for the container to identify the directory path and for snakemake to see the output files. 

```
cd /home/tiakju/test_folder/steinbock_snakemake # path to where steinbock-snakemake is cloned

snakemake -c 4 --use-singularity --singularity-args "-B /home/tiakju/mnt/crucial/SynologyDrive/pdac_projects/steinbock_snakemake" --configfile data/test_mcd/test_mcd.yaml
```

### Inspect pipeline outputs

This concludes this tutorial!

```
└── 📁test_mcd
    └── 📁deepcell # outputs for deepcell segmentation and masks per ROI
        └── 📁intensities
            └── test_018.csv
        └── 📁neighbors
            └── test_018.csv
        └── 📁nuclei
            └── test_018.tiff
        └── 📁regionprops
            └── test_018.csv
        └── 📁whole_cell
            └── test_018.tiff
    └── 📁export # Exports per channel tiff for each ROI, an ome.tiff file and a anndata object
        └── 📁test_018
            └── ArAr80_80ArAr.tiff
            └── Dy162_162Dy_h5454_Chr10SAT.tiff
            └── Er166_166Er_h3838_Chr1SAT.tiff
            └── Er167_167Er_h3838_Chr1SAT.tiff
            └── Gd158_158Gd_h5050_Chr2SAT.tiff
            └── Ir191_191Ir_DNA1.tiff
            └── Ir193_193Ir_DNA2.tiff
            └── Pb206_206Pb.tiff
            └── Pb208_208Pb.tiff
            └── Xe126_126Xe.tiff
            └── Xe131_131Xe.tiff
            └── Xe134_134Xe.tiff
            └── test_018_mask.tiff
        └── test_018.ome.tiff
        └── test_mcd.h5ad
    └── 📁img # Multichannel Tiff folder
        └── images.csv
        └── test_018.tiff
    └── 📁mcd
        └── test.mcd
    └── panel.csv # panel file containing all markers
    └── panel_deepcell.csv # panel file containing markers for nuclear and cytoplasm segmentation
    └── test_mcd.yaml
```