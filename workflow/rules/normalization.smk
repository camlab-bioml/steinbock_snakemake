projects, ROI = glob_wildcards(rules.generate_tiff.output.place + "/" + "{ROI}" + ".tiff")

normalization_output = {
    "tiffs_arcsinh": expand("data/{projects}/img/arcsinh/{ROI}.tiff", projects = projects, ROI = ROI),
    "tiffs_zscale": expand("data/{projects}/img/zscale/{ROI}.tiff", projects = projects, ROI = ROI),
    "plots": expand("data/{projects}/img/plots/{ROI}", projects = projects, ROI = ROI)
}

rule arcsinh_normalization:
    input:
        raw = "data/{projects}/img/raw/{ROI}.tiff"
    output:
        normalized = "data/{projects}/img/arcsinh/{ROI}.tiff"
    threads: 1
    conda: "steinbock-snakemake"
    script:
        "normalization/arcsinh_normalization.py"

rule zscale_normalization:
    input:
        raw = "data/{projects}/img/raw/{ROI}.tiff"
    output:
        normalized = "data/{projects}/img/zscale/{ROI}.tiff"
    params:
        alpha = 1,
        beta = 0,
    threads: 1
    conda: "steinbock-snakemake"
    script:
        "normalization/zscale_normalization.py"

rule plot_normalization:
    input:
        raw = "data/{projects}/img/raw/{ROI}.tiff",
        arcsinh = rules.arcsinh_normalization.output.normalized,
        zscale = rules.zscale_normalization.output.normalized,
        panel = rules.create_panel.output.panel
    output:
        plots = directory("data/{projects}/img/plots/{ROI}")
    threads: 1
    conda: "steinbock-snakemake"
    script:
        "normalization/plot_normalization.py"
