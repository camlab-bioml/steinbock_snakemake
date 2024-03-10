rule cellprofiler_prepare:
    input:
        rules.ilastik_prepare.output
    output:
        "data/{projects}/cellprofiler/cell_segmentation.cppipe"
    singularity:
        "workflow/envs/steinbock-gpu.sif"
    shell:
        """
            steinbock segment cellprofiler prepare \
                -o {output} 
        """