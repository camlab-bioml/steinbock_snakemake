preprocessing_output = {
    "panel": expand("data/{projects}/panel.csv", projects = projects),
    "tiff_folder": expand("data/{projects}/img/raw", projects = projects),
    "tiff_metadata": expand("data/{projects}/img/images.csv", projects = projects)
}

rule create_panel:
    input:
        a = "data/{projects}/mcd"
    output:
        "data/{projects}/panel.csv"
    conda: "steinbock-snakemake"
    shell:
        """
        steinbock preprocess imc panel \
            --mcd {input.a} \
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
        place = directory("data/{projects}/img/raw"),
        stats = "data/{projects}/img/images.csv"
    conda: "steinbock-snakemake"
    shell:
        """
        steinbock preprocess imc images --mcd {input.b} \
                                        --txt {input.b} \
                                        --panel {input.a} \
                                        --imgout {output.place} \
                                        --hpf {params.hpf} \
                                        --infoout {output.stats} \
                                        --verbosity DEBUG \
                                        --strict True
        """
