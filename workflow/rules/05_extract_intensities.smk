rule extract_intensities:
    input:
        masks = rules.deepcell_wholecell.output,
        img = rules.generate_tiff.output.place,
        panel = rules.create_panel.output
    params:
        aggr = config["aggr"]
    output:
        directory("data/{projects}/deepcell/intensities")
    singularity:
        config["container"]
    threads: 24
    shell:
        """
        steinbock measure intensities \
            --img {input.img} \
            --masks {input.masks} \
            --panel {input.panel} \
            --aggr {params.aggr} \
            -o {output} \
            -v DEBUG
        """