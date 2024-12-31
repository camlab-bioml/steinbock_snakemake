export = {
    "export_all_zarr": expand("data/{projects}/export/{projects}.zarr", projects = projects),
    "export_all_h5ad": expand("data/{projects}/export/{projects}.h5ad", projects = projects),
    "export_all_ome": expand("data/{projects}/export/ome", projects = projects),
    "export_all_umap": expand("data/{projects}/export/umap/umap_min_dist_{umap_dist}_coordinates.csv", projects = projects, umap_dist = umap_dist),
    "export_umap_plots": expand("data/{projects}/export/umap/umap_min_dist_{umap_dist}.png", projects = projects, umap_dist = umap_dist),

}

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
        h5ad = "data/{projects}/export/{projects}_export.h5ad",
        zarr = directory("data/{projects}/export/{projects}.zarr")
    threads: 24
    shell:
        """
        steinbock export ome --img {input.img} --panel {input.panel} -o {output.dir}
        steinbock export histocat --img {input.img} --masks {input.masks} --panel {input.panel} -o {output.dir}
        steinbock export anndata --intensities {input.intensities} --data {input.regionprops} --neighbors {input.neighbors} --panel {input.panel} --info {input.imginfo} --format zarr -o {output.zarr}
        steinbock export anndata --intensities {input.intensities} --data {input.regionprops} --neighbors {input.neighbors} --panel {input.panel} --info {input.imginfo} --format h5ad -o {output.h5ad}
        """

rule phenograph:
    input:
        h5ad = rules.export_all.output.h5ad
    params:
        script = srcdir("export/phenograph_clustering.py"),
        k = config["phenograph_k"],
        min_cluster_size = config["phenograph_min_cluster_size"]
    output:
        h5ad = "data/{projects}/export/{projects}.h5ad"
    shell:
        """
        mv {input.h5ad} {output.h5ad} && python {params.script} --input {output.h5ad} --output {output.h5ad} -k {params.k} -m {params.min_cluster_size} --normalize
        """

rule umap:
    input:
        h5ad = rules.phenograph.output.h5ad
    params:
        script = srcdir("export/dimension_reduction.py"),
        min_dist = lambda wildcards: wildcards.umap_dist
    conda: "steinbock-snakemake"
    output:
        coords = "data/{projects}/export/umap/umap_min_dist_{umap_dist}_coordinates.csv",
        plots = "data/{projects}/export/umap/umap_min_dist_{umap_dist}.png"
    shell:
        """
        mkdir -p data/{projects}/export/umap/
        cp {input.h5ad} data/{projects}/export/tmp_{params.min_dist}.h5ad
        python {params.script} --input data/{projects}/export/tmp_{params.min_dist}.h5ad --output-coords {output.coords} \
        --normalize --min-dist {params.min_dist} --plot --output-plot {output.plots}
        rm data/{projects}/export/tmp_{params.min_dist}.h5ad
        """