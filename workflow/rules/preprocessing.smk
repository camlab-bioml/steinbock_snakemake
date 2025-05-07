version: "0.0.5"

preprocessing_output = {
    "panel": expand("data/{projects}/panel.csv", projects = projects),
    "tiff_folder": expand("data/{projects}/img/raw", projects = projects),
    "tiff_metadata": expand("data/{projects}/img/images.csv", projects = projects)
}

rule create_panel:
    input:
        mcd = "data/{projects}/mcd"
    output:
        panel = "data/{projects}/panel.csv"
    conda: "steinbock-snakemake"
    shell:
        """
        steinbock preprocess imc panel \
            --mcd {input.mcd} \
            -o {output} \
            --verbosity INFO
        """

rule generate_tiff:
    input:
        panel = rules.create_panel.output,
        mcd = "data/{projects}/mcd"
    params:
        hpf = config['hpf']
    output:
        tiff_folder = directory("data/{projects}/img/raw"),
        tiff_metadata = "data/{projects}/img/images.csv"
    conda: "steinbock-snakemake"
    shell:
        """
        steinbock preprocess imc images --mcd {input.mcd} \
                                        --txt {input.mcd} \
                                        --panel {input.panel} \
                                        --imgout {output.tiff_folder} \
                                        --hpf {params.hpf} \
                                        --infoout {output.tiff_metadata} \
                                        --verbosity DEBUG \
                                        --strict True
        """
