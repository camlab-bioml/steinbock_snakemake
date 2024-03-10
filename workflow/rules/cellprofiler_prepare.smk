rule cellprofiler_prepare:
    input:
        rules.ilastik.output
    output:
        "data/{projects}/cellprofiler/cell_segmentation.cppipe"
    shell:
        """
            steinbock segment cellprofiler prepare \
                -o {output} 
        """