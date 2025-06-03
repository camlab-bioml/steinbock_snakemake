deepcell = {
    "panel_deepcell": expand("data/{projects}/panel_deepcell.csv", projects = projects),
    "deepcell_nuclei": expand("data/{projects}/deepcell/nuclei", projects = projects),
    "deepcell_wholecell": expand("data/{projects}/deepcell/cell", projects = projects),
    "deepcell_overlay": expand("data/{projects}/deepcell/overlay", projects = projects)
}

rule deepcell_prepare:
    input:
        rules.create_panel.output
    output:
        panel_deepcell = "data/{projects}/panel_deepcell.csv"
    params:
        script = "workflow/scripts/deepcell_prepare.py",
        cytoplasm = config["cytoplasm"],
        nuclear = config["nuclear"]
    conda: "steinbock-snakemake"
    shell:
        """
        python {params.script} \
                --panel {input} \
                --nuclear {params.nuclear} \
                --cytoplasm {params.cytoplasm} \
                --output {output}
        """

rule deepcell_nuclei:
    input:
        img = rules.generate_tiff.output.tiff_folder,
        panel = rules.deepcell_prepare.output.panel_deepcell
    params:
        app = config["deepcell_app"],
        model = config["deepcell_model"],
        model_path = config["deepcell_modelpath"],
        px_size = config["deepcell_pxsize"]
    output:
        directory("data/{projects}/deepcell/nuclei")
    conda: "steinbock-snakemake"
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
            --panel {input.panel} \
            --pixelsize {params.px_size} \
            -o {output} \
            -v DEBUG
        """

rule deepcell_wholecell:
    input:
        img = rules.generate_tiff.output.tiff_folder,
        panel = rules.deepcell_prepare.output
    params:
        app = config["deepcell_app"],
        model = config["deepcell_model"],
        model_path = config["deepcell_modelpath"],
        px_size = config["deepcell_pxsize"]
    output:
        directory("data/{projects}/deepcell/cell")
    conda: "steinbock-snakemake"
    threads: 24
    shell:
        """
        steinbock segment deepcell \
            --app {params.app} \
            --model {params.model} \
            --modeldir {params.model_path} \
            --type whole-cell \
            --img {input.img} \
            --minmax \
            --panel {input.panel} \
            --pixelsize {params.px_size} \
            -o {output} \
            -v DEBUG
        """

rule nuclear_overlay:
    input:
        img = rules.generate_tiff.output.tiff_folder,
        panel = rules.deepcell_prepare.output.panel_deepcell,
        mask_dir = rules.deepcell_wholecell.output
    output:
        directory("data/{projects}/deepcell/overlay")
    params:
        script = "workflow/scripts/nuclear_overlay.py"
    shell:
        """
        mkdir -p {output}
        python {params.script} -i {input.img} -p {input.panel} -o {output} -m {input.mask_dir} -t tiff
        """
