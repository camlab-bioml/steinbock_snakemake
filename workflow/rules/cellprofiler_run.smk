rule cellprofiler_run:
    input:
        prepare = rules.cellprofiler_prepare.output,
        probabs = rules.ilastik_run.output.probabilities
    params:
        python_path = "",
        cellprofiler_module = "",
        plugins_dir = ""
    output:
        "data/{projects}/cellprofiler/masks"
    shell:
        """
            steinbock segment cellprofiler run \
                    --python \
                    --cellprofiler \
                    --plugins-directory \
                    --pipe {input.prepare} \
                    --probabs {input.probabs} \
                    -o {output}
        """