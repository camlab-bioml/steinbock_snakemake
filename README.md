# Steinbock-snakemake v0.0.1

## Set-up pipeline environment (via conda, pipenv)

We will first set up the pipeline environment via conda. This will install all dependencies including snakemake and singularity for use.

```
git clone <steinbock-snakemake> && cd steinbock-snakemake # clone git repo, then move into the pipeline folder. 
conda create env -y ---name steinbock-snakemake --file=workflow/env/environment.yml
conda activate steinbock-snakemake
```

Alternatively we can use pipenv to just install snakemake given support for singularity in a cluster environment. (Not recoomended)
```
pipenv install -r workflow/env/requirements.txt 
```

## Clone docker containers

### Cloning CPU-based container (Recommended)
After successful creation of the environment, we will then generate a singularity container for use with the snakemake pipeline. This container will require about ~4GB of storage, and will be located under workflow/envs

```
singularity pull docker://ghcr.io/bodenmillergroup/steinbock:0.16.1 workflow/envs/steinbock-cpu.sif
```

### Cloning GPU-based container (Not Recommended)
This will clone the GPU container for steinbock. Note, you will require `nvidia-container-toolkit` and require CUDA drivers to run the pipeline. Only use this for large datasets that require GPU support.
```
singularity pull docker://ghcr.io/bodenmillergroup/steinbock:0.16.1-gpu workflow/envs/steinbock-cpu.sif
```

## Running Steinbock-snakemake

### Data input and output structures
Now that the containers and environment is set-up, we will run and test the pipeline on a small mcd file. The structure of a completed pipeline:
```
── 📁data
    └── 📁test_mcd
        └── 📁deepcell # Folder for all deepcell outputs
        └── 📁export # Folder for steinbock outputs. Will output a h5ad per ROI.
        └── 📁img # Folder for tiff image stacks
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
cytoplasm: "Gd160 Ho165"
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
    - myproject # add the folder name to your project
```

### Channels for segmenting cells
To help mesmer segment cells, steinbock will require the user to input the mass channels for cytoplasm specific markers (e.g panCK, Actin etc). You should include those channels under `cytoplasm` in Rule Group 4.
```
cytoplasm: "Gd160 Ho165" # <-- add mass channels that corresspond to the cytoplasm
deepcell_app: "mesmer"
deepcell_model: "MultiplexSegmentation"
deepcell_modelpath: "/opt/keras/models"
deepcell_pxsize: 1
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
After editing `config.yml`, we would want to ensure we use the settings for each project. A simple solution is to copy the `config.yml` file over to the project directory before calling the snakemake pipeline. 
```
cp config/config.yaml data/test_mcd/test_mcd.yaml
```

### Running the snakemake pipeline
After completing the preflight configuration, we can run the pipeline. The `snakemake` command is invoked in the `steinbock-snakemake` directory:
1. Specify cores required using the `-c` flag.
2. Specify the configuration file used using the `--configfile` flag.
3. Specify bind mount for singularity container using the `--singularity-args` flag. This is required for the container to identify the directory path and for snakemake to see the output files. 

```
cd ~/steinbock-snakemake/ # path to where steinbock-snakemake is cloned

snakemake -c 4 \
          --use-singularity \
          --singularity-args "-B ~/steinbock_snakemake" \
          --configfile data/test_mcd/test_mcd.yaml \
```

