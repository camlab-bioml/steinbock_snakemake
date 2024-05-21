mesmer_pipeline_outputs = {
    panel: expand("data/{projects}/panel.csv", projects = projects),
    tiff_folder: directory(expand("data/{projects}/img", projects = projects)),
    tiff_metadata: expand("data/{projects}/img/images.csv", projects = projects),
    panel_deepcell: expand("data/{projects}/panel_deepcell.csv", projects = projects),
    deepcell_nuclei: directory(expand("data/{projects}/deepcell/nuclei"), projects = projects),
    deepcell_wholecell: directory(expand("data/{projects}/deepcell/whole_cell", projects = projects)),
    extract_intensities: directory(expand("data/{projects}/deepcell/intensities", projects = projects)),
    extract_neighbors: directory(expand("data/{projects}/deepcell/neighbors", projects = projects)),
    extract_regionprops: directory(expand("data/{projects}/deepcell/regionprops", projects = projects)),
    export_all_zarr: directory(expand("data/{projects}/export/{projects}.zarr"), projects = projects),
    export_all_h5ad: expand("data/{projects}/export/{projects}.h5ad", projects = projects),
    export_all_ome: directory(expand("data/{projects}/export/ome"), projects = projects)
}

rule create_panel:
    input:
        "data/{projects}/mcd"
    output:
        "data/{projects}/panel.csv"
    singularity:
        config["container"]
    shell:
        """
        steinbock preprocess imc panel \
            --mcd {input} \
            -o {output} \
            --verbosity INFO
        """

rule generate_tiff:
    input:
        a = rules.create_panel.output,
        b = "data/{projects}/mcd"
    params:
        hpf = config['hpf']
    output:
        place = directory("data/{projects}/img"),
        stats = "data/{projects}/img/images.csv"
    singularity:
        config["container"]
    shell:
        """
        steinbock preprocess imc images --mcd {input.b} \
                                        --txt {input.b} \
                                        --panel {input.a} \
                                        --imgout {output.place} \
                                        --hpf {params.hpf} \
                                        --infoout {output.stats} \
                                        --verbosity DEBUG
        """

rule deepcell_prepare:
    input:
        "data/{projects}/panel.csv"
    output:
        "data/{projects}/panel_deepcell.csv"
    params:
        script = "workflow/scripts/deepcell_prepare.py",
        cytoplasm = config["cytoplasm"],
        nuclei = config["nuclear"]
    singularity:
        config["container"]
    shell:
        """
        python {params.script} \
                --workdir data/{projects}/ \
                --nuclear {params.cytoplasm} \
                --cytoplasm {params.cytoplasm} \
                --output panel_deepcell.csv
        """

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
        config["container"]
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

rule deepcell_wholecell:
    input:
        img = rules.generate_tiff.output.place,
        panel = rules.deepcell_prepare.output
    params:
        app = config["deepcell_app"],
        model = config["deepcell_model"],
        model_path = config["deepcell_modelpath"],
        px_size = config["deepcell_pxsize"],
    output:
        directory("data/{projects}/deepcell/whole_cell")
    singularity:
        config["container"]
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
            --zscore \
            --panel {input.panel} \
            --pixelsize {params.px_size} \
            -o {output} \
            -v DEBUG
        """

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

rule export_all:
    input:
        masks = rules.deepcell_wholecell.output,
        img = rules.generate_tiff.output.place,
        imginfo = rules.generate_tiff.output.stats,
        panel = rules.create_panel.output,
        intensities = rules.extract_intensities.output,
        regionprops = rules.extract_regionprops.output,
        neighbors = rules.extract_neighbors.output
    params:
        aggr = config["aggr"]
    output:
        dir = directory("data/{projects}/export/ome"),
        h5ad = "data/{projects}/export/{projects}.h5ad",
        zarr = directory("data/{projects}/export/{projects}.zarr")
    singularity:
        config["container"]
    threads: 24
    shell:
        """
        steinbock export ome --img {input.img} --panel {input.panel} -o {output.dir}
        steinbock export histocat --img {input.img} --masks {input.masks} --panel {input.panel} -o {output.dir}
        steinbock export anndata --intensities {input.intensities} --data {input.regionprops} --neighbors {input.neighbors} --panel {input.panel} --info {input.imginfo} --format zarr -o {output.zarr}
        steinbock export anndata --intensities {input.intensities} --data {input.regionprops} --neighbors {input.neighbors} --panel {input.panel} --info {input.imginfo} --format h5ad -o {output.h5ad}
        """