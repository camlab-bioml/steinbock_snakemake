rule create_panel:
    input:
        "data/{projects}/mcd"
    output:
        "data/{projects}/panel.csv"
    singularity:
        "workflow/envs/steinbock-gpu.sif"
    shell:
        """
        steinbock preprocess imc panel \
            --mcd {input} \
            -o {output} \
            --verbosity INFO
        """