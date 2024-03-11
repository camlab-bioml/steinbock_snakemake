rule extract_neighbors:
    input:
        masks = rules.deepcell_wholecell.output
    params:
        type = config["neighbor_type"],
        dmax = config["dmax"],
        kmax = config["kmax"]
    output:
        directory("data/{projects}/deepcell/neighbors")
    singularity:
        config["container"]
    threads: 24
    shell:
        """
        steinbock measure neighbors \
            --masks {input.masks} \
            --type {params.type} \
            --dmax {params.dmax} \
            --kmax {params.kmax} \
            -o {output} \
            -v DEBUG
        """