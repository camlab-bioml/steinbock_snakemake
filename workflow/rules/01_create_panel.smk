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