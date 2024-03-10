rule deepcell_wholecell:
    input:
        img = rules.generate_tiff.output.place,
        panel = rules.create_panel.output
    params:
        app = config["deepcell_app"],
        model = config["deepcell_model"],
        model_path = config["deepcell_modelpath"],
        px_size = config["deepcell_pxsize"]
    output:
        "data/{projects}/deepcell/whole_cell"
    shell:
        """
        steinbock segment deepcell \
            --app {params.app} \
            --model {params.model} \
            --type whole-cell \
            --img {input.img} \
            --minmax \
            --zscore \
            --panel {input.panel} \
            --pixelsize {params.px_size}
            --preprocess \
            --postprocess \
            -o {output} \
            -v DEBUG
        """