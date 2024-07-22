normalization_output = {
     "arcsinh": expand("data/{projects}/img/arcsinh/", projects = projects),
     "minmax": expand("data/{projects}/img/minmax/", projects = projects),
     "zscale": expand("data/{projects}/img/zscale/", projects = projects),
   #  "preview": expand("data/{projects}/img/preview", projects = projects)
}

rule arcsinh:
    input:
        raw_folder = rules.generate_tiff.output.place
    output:
        arcsinh_folder = directory("data/{projects}/img/arcsinh/")
    threads: 1
    conda: "steinbock-snakemake"
    script:
        "normalization/arcsinh_normalization.py"

rule minmax:
    input:
        raw_folder = rules.generate_tiff.output.place
    output:
        minmax_folder = directory("data/{projects}/img/minmax/")
    threads: 1
    conda: "steinbock-snakemake"
    script:
        "normalization/minmax_normalization.py"

rule zscale:
    input:
        raw_folder = rules.generate_tiff.output.place
    output:
        zscale_folder = directory("data/{projects}/img/zscale/")
    threads: 1
    conda: "steinbock-snakemake"
    script:
        "normalization/zscale_normalization.py"

# rule preview:
#     input:
#         raw = rules.generate_tiff.output.place,
#         arcsinh = rules.arcsinh.output.arcsinh_folder,
#         zscale = rules.zscale.output.zscale_folder,
#         minmax = rules.minmax.output.minmax_folder, 
#         panel = rules.create_panel.output,
#         images = rules.generate_tiff.output.stats
#     output:
#         directory("data/{projects}/img/preview")
#     threads: 1
#     conda: "steinbock-snakemake"
#     script:
#         "normalization/preview.py"