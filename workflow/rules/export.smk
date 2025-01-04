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
        h5ad = temp("{projects}_export.h5ad"),
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
        h5ad = rules.export_all.output.h5ad,
        panel = rules.create_panel.output
    params:
        script = srcdir("export/phenograph_clustering.py"),
        k = config["phenograph_k"],
        min_cluster_size = config["phenograph_min_cluster_size"],
        channels_ignore = channels_ignore
    output:
        h5ad = temp("{projects}_pheno.h5ad")
    shell:
        """
        python {params.script} -i {input.h5ad} -o {output.h5ad} -k {params.k} -m {params.min_cluster_size} --normalize -ci {params.channels_ignore} -p {input.panel}
        """

rule umap:
    input:
        h5ad = rules.phenograph.output.h5ad,
        panel = rules.create_panel.output
    params:
        script = srcdir("export/dimension_reduction.py"),
        min_dist = lambda wildcards: wildcards.umap_dist,
        tmp_file = temp("tmp_{params.min_dist}.h5ad"),
        channels_ignore = channels_ignore
    output:
        coords = "data/{projects}/export/umap/umap_min_dist_{umap_dist}_coordinates.csv",
        plots = "data/{projects}/export/umap/umap_min_dist_{umap_dist}.png",
    shell:
        """
        mkdir -p data/{projects}/export/umap/
        cp {input.h5ad} {params.tmp_file}
        python {params.script} --i {params.tmp_file} -oc {output.coords} \
        --normalize -md {params.min_dist} --plot -op {output.plots} -ci {params.channels_ignore} -p {input.panel}
        """
        
rule collect_umap_dists:
    input:
        export_h5ad = rules.export_all.output.h5ad,
        pheno_h5ad = rules.phenograph.output.h5ad,
        plots = expand(
            "data/{projects}/export/umap/umap_min_dist_{umap_dist}.png",
            projects=projects, umap_dist=umap_dist)
    output:
        h5ad = "data/{projects}/export/{projects}.h5ad"
    params:
        script = srcdir("export/collect_umap_dists.py"),
        dir = "data/{projects}/export/umap/"
    shell:
        """
        python {params.script} -i {input.pheno_h5ad} -cd {params.dir} -o {output.h5ad}
        """