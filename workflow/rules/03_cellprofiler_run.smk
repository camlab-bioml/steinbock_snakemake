rule cellprofiler_run:
    input:
        prepare = rules.cellprofiler_prepare.output,
        probabs = rules.ilastik_run.output.probabilities
    params:
        python_path = config["python_path"],
        plugins_dir = config["plugins_path"]
    output:
        "~/steinbock_snakemake/data/{projects}/cellprofiler/masks"
    singularity:
        config["container"]
    shell:
        """
            steinbock segment cellprofiler run \
                    --python {params.python_path} \
                    --plugins-directory {params.plugins_dir} \
                    --pipe {input.prepare} \
                    --probabs {input.probabs} \
                    -o {output}
        """