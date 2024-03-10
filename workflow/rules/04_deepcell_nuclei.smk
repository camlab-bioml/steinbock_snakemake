rule deepcell_nuclei:
    input:
        img = rules.generate_tiff.output.place,
        panel = rules.deepcell_prepare.output
    params:
        app = config["deepcell_app"],
        model = config["deepcell_model"],
        model_path = config["deepcell_modelpath"],
        px_size = config["deepcell_pxsize"]
    output:
        directory("data/{projects}/deepcell/nuclei")
    singularity:
        "workflow/envs/steinbock-gpu.sif"
    threads: 24
    shell:
        """
        steinbock segment deepcell \
            --app {params.app} \
            --model {params.model} \
            --modeldir {params.model_path} \
            --type nuclear \
            --img {input.img} \
            --minmax \
            --zscore \
            --panel {input.panel} \
            --pixelsize {params.px_size} \
            -o {output} \
            -v DEBUG
        """