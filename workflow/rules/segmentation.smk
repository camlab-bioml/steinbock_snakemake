deepcell = {
    "panel_deepcell": expand("data/{projects}/panel_deepcell.csv", projects = projects),
    "deepcell_nuclei": expand("data/{projects}/deepcell/nuclei", projects = projects),
    "deepcell_wholecell": expand("data/{projects}/deepcell/cell", projects = projects),
    "deepcell_overlay": expand("data/{projects}/deepcell/overlay", projects = projects),
    "match_segmentation_files_for_quantification": expand("data/{projects}/logs/clean_unmatched_files.done", projects = projects)
}

rule deepcell_prepare:
    input:
        rules.create_panel.output if not process_tiff else config['panel_tiff']
    output:
        panel_deepcell = "data/{projects}/panel_deepcell.csv"
    params:
        script = "workflow/scripts/deepcell_prepare.py",
        cytoplasm = config["cytoplasm"],
        nuclear = config["nuclear"]
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
        img = rules.generate_tiff.output.tiff_folder if not process_tiff else config['tiff_path'],
        panel = rules.deepcell_prepare.output.panel_deepcell
    params:
        app = config["deepcell_app"],
        model = config["deepcell_model"],
        model_path = config["deepcell_modelpath"],
        px_size = config["deepcell_pxsize"],
        min_max = "--minmax" if not process_tiff else ""
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
            {params.min_max} \
            --panel {input.panel} \
            --pixelsize {params.px_size} \
            -o {output} \
            -v DEBUG
        """

rule deepcell_wholecell:
    input:
        img = rules.generate_tiff.output.tiff_folder if not process_tiff else config['tiff_path'],
        panel = rules.deepcell_prepare.output.panel_deepcell
    params:
        app = config["deepcell_app"],
        model = config["deepcell_model"],
        model_path = config["deepcell_modelpath"],
        px_size = config["deepcell_pxsize"],
        min_max = "--minmax" if not process_tiff else ""
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
            {params.min_max} \
            --panel {input.panel} \
            --pixelsize {params.px_size} \
            -o {output} \
            -v DEBUG
        """

rule match_segmentation_files_for_quantification:
    input:
        dir = rules.deepcell_wholecell.output
    output:
        out_log = temp("data/{projects}/logs/clean_unmatched_files.done")
    params:
        raw_tiff = rules.generate_tiff.output.tiff_folder if not process_tiff else config['tiff_path']
    conda: "steinbock-snakemake"
    threads: 2
    shell:
        """
        mkdir -p "$(dirname {params.raw_tiff})/raw_not_quantified" && \
        comm -23 <(ls {params.raw_tiff} | sort) <(ls {input.dir} | sort) \
        | xargs -r -I{{}} mv "{params.raw_tiff}/{{}}" "$(dirname {params.raw_tiff})/raw_not_quantified/"
        touch {output}
        """

rule nuclear_overlay:
    input:
        img = rules.generate_tiff.output.tiff_folder if not process_tiff else config['tiff_path'],
        panel = rules.deepcell_prepare.output.panel_deepcell,
        mask_dir = rules.deepcell_wholecell.output
    output:
        directory("data/{projects}/deepcell/overlay")
    params:
        script = "workflow/scripts/nuclear_overlay.py"
    conda: "steinbock-snakemake"
    threads: 2
    shell:
        """
        mkdir -p {output}
        python {params.script} -i {input.img} -p {input.panel} -o {output} -m {input.mask_dir} -t tiff
        """
