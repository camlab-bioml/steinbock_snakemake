rule cellprofiler_prepare:
    input:
        rules.ilastik_prepare.output
    output:
        "data/{projects}/cellprofiler/cell_segmentation.cppipe"
    singularity:
        config["container"]
    shell:
        """
            steinbock segment cellprofiler prepare \
                -o {output} 
        """