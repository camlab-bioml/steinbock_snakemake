export = {
    "export_all_zarr": expand("data/{projects}/export/{projects}.zarr", projects = projects),
    "export_all_h5ad": expand("data/{projects}/export/{projects}.h5ad", projects = projects),
    "export_all_ome": expand("data/{projects}/export/ome", projects = projects)
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
        h5ad = "data/{projects}/export/{projects}.h5ad",
        zarr = directory("data/{projects}/export/{projects}.zarr")
    threads: 24
    shell:
        """
        steinbock export ome --img {input.img} --panel {input.panel} -o {output.dir}
        steinbock export histocat --img {input.img} --masks {input.masks} --panel {input.panel} -o {output.dir}
        steinbock export anndata --intensities {input.intensities} --data {input.regionprops} --neighbors {input.neighbors} --panel {input.panel} --info {input.imginfo} --format zarr -o {output.zarr}
        steinbock export anndata --intensities {input.intensities} --data {input.regionprops} --neighbors {input.neighbors} --panel {input.panel} --info {input.imginfo} --format h5ad -o {output.h5ad}
        """