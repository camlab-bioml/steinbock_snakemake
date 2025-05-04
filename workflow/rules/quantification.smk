quantification = {
    "extract_intensities": expand("data/{projects}/quantification/intensities", projects = projects),
    "extract_neighbors": expand("data/{projects}/quantification/neighbors", projects = projects),
    "extract_regionprops": expand("data/{projects}/quantification/regionprops", projects = projects),
}

rule extract_intensities:
    input:
        masks = rules.deepcell_wholecell.output,
        img = rules.generate_tiff.output.tiff_folder,
        panel = rules.create_panel.output.panel
    params:
        aggr = config["aggr"]
    output:
        directory("data/{projects}/quantification/intensities")
    threads: 24
    conda: 'steinbock-snakemake'
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

rule extract_neighbors:
    input:
        masks = rules.deepcell_wholecell.output
    params:
        type = config["neighbor_type"],
        dmax = config["dmax"],
        kmax = config["kmax"]
    output:
        directory("data/{projects}/quantification/neighbors")
    threads: 24
    conda: 'steinbock-snakemake'
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

rule extract_regionprops:
    input:
        masks = rules.deepcell_wholecell.output,
        img = rules.generate_tiff.output.tiff_folder
    params:
        aggr = config["aggr"]
    output:
        directory("data/{projects}/quantification/regionprops")
    threads: 24
    conda: 'steinbock-snakemake'
    shell:
        """
        steinbock measure regionprops \
            --img {input.img} \
            --masks {input.masks} \
            -o {output} \
            -v DEBUG
        """