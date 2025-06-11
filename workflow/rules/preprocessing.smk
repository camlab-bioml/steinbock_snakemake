version: "0.0.6"

preprocessing_output = {
    "panel": expand("data/{projects}/panel.csv", projects = projects) if not process_tiff else [],
    "tiff_folder": expand("data/{projects}/img/raw", projects = projects) if not process_tiff else [],
    "tiff_metadata": expand("data/{projects}/img/images.csv", projects = projects) if not process_tiff else []
}

rule create_panel:
    input:
        mcd = "data/{projects}/mcd"
    output:
        panel = "data/{projects}/panel.csv" if not process_tiff else []
    conda:
        "steinbock-snakemake"
    params:
        type = "mcd" if not process_txt else "txt"
    shell:
        """
        steinbock preprocess imc panel \
                --{params.type} {input.mcd} \
                -o {output.panel} \
                --verbosity INFO
        """

rule generate_tiff:
    input:
        panel = rules.create_panel.output,
        mcd = "data/{projects}/mcd" if not process_tiff else []
    params:
        hpf = config['hpf'],
        skip = process_tiff
    output:
        tiff_folder = directory("data/{projects}/img/raw") if not process_tiff else [],
        tiff_metadata = "data/{projects}/img/images.csv" if not process_tiff else []
    conda:
        "steinbock-snakemake"
    shell:
        """
        steinbock preprocess imc images \
                --mcd {input.mcd} \
                --txt {input.mcd} \
                --panel {input.panel} \
                --imgout {output.tiff_folder} \
                --hpf {params.hpf} \
                --infoout {output.tiff_metadata} \
                --verbosity DEBUG \
                --strict True
        """
