**TO DO ADD INSTRUCTIONS

To run the pipeline, generate a project folder under data/{your_project} and place raw mcd files in data/{your_project}/mcd.

The pipeline uses mesmer and requires input on cytoplasm markers. To specify them, edit the "cytoplasm" config in config/config.yaml. 

To run the container, use the following command
snakemake -c 24 --use-singularity --singularity-args "-B <your local mount>"

The singularity containser is found under workflow/envs/steinbock-gpu.sif

The output will be populated to data/{your_project}.

General Pipeline
extract panel -> filter hot pixels -> mesmer segmentation -> extract single cell intensities -> measure regions -> find neighbors -> export anndata
