quantification = {
    "extract_intensities": expand("data/{projects}/quantification/intensities", projects = projects),
    "extract_neighbors": expand("data/{projects}/quantification/neighbors", projects = projects),
    "extract_regionprops": expand("data/{projects}/quantification/regionprops", projects = projects),
    "intensities_nuclei": expand("data/{projects}/quantification/nuclei/intensities", projects = projects) if quantify_nuclei else [],
    "neighbors_nuclei": expand("data/{projects}/quantification/nuclei/neighbors", projects = projects) if quantify_nuclei else [],
    "regionprops_nuclei": expand("data/{projects}/quantification/nuclei/regionprops", projects = projects) if quantify_nuclei else []
}

rule extract_intensities:
    input:
        masks = rules.deepcell_wholecell.output,
        img = rules.generate_tiff.output.tiff_folder if not process_tiff else config['tiff_path'],
        panel = rules.deepcell_prepare.output.panel_deepcell if not process_tiff else config['panel_tiff'],
        log = rules.match_segmentation_files_for_quantification.output.out_log
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
        img = rules.generate_tiff.output.tiff_folder if not process_tiff else config['tiff_path'],
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

rule extract_intensities_nuclei:
    input:
        masks = rules.deepcell_nuclei.output,
        img = rules.generate_tiff.output.tiff_folder if not process_tiff else config['tiff_path'],
        panel = rules.deepcell_prepare.output.panel_deepcell if not process_tiff else config['panel_tiff'],
        log = rules.match_segmentation_files_for_quantification.output.out_log
    params:
        aggr = config["aggr"],
        quantify = quantify_nuclei
    output:
        directory("data/{projects}/quantification/nuclei/intensities")
    threads: 24
    conda: "steinbock-snakemake"
    shell:
        r"""
        if [ "{params.quantify}" = "True" ]; then
            steinbock measure intensities \
                --img {input.img} \
                --masks {input.masks} \
                --panel {input.panel} \
                --aggr {params.aggr} \
                -o {output} \
                -v DEBUG
        fi
        """


rule extract_neighbors_nuclei:
    input:
        masks = rules.deepcell_nuclei.output
    params:
        type = config["neighbor_type"],
        dmax = config["dmax"],
        kmax = config["kmax"],
        quantify = quantify_nuclei
    output:
        directory("data/{projects}/quantification/nuclei/neighbors")
    threads: 24
    conda: 'steinbock-snakemake'
    shell:
        r"""
        if [ "{params.quantify}" = "True" ]; then
            steinbock measure neighbors \
            --masks {input.masks} \
            --type {params.type} \
            --dmax {params.dmax} \
            --kmax {params.kmax} \
            -o {output} \
            -v DEBUG
        fi
        """

rule extract_regionprops_nuclei:
    input:
        masks = rules.deepcell_nuclei.output,
        img = rules.generate_tiff.output.tiff_folder if not process_tiff else config['tiff_path'],
    params:
        aggr = config["aggr"],
        quantify = quantify_nuclei
    output:
        directory("data/{projects}/quantification/nuclei/regionprops")
    threads: 24
    conda: 'steinbock-snakemake'
    shell:
        r"""
        if [ "{params.quantify}" = "True" ]; then
            steinbock measure regionprops \
            --img {input.img} \
            --masks {input.masks} \
            -o {output} \
            -v DEBUG
        fi
        """