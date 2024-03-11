rule extract_regionprops:
    input:
        masks = rules.deepcell_wholecell.output,
        img = rules.generate_tiff.output.place
    params:
        aggr = config["aggr"]
    output:
        directory("data/{projects}/deepcell/regionprops")
    singularity:
        config["container"]
    threads: 24
    shell:
        """
        steinbock measure regionprops \
            --img {input.img} \
            --masks {input.masks} \
            -o {output} \
            -v DEBUG
        """