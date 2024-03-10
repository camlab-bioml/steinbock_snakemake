rule create_panel:
    input:
        "data/{projects}/mcd"
    output:
        "data/{projects}/panel.csv"
    shell:
        """
        steinbock preprocess imc panel \
            -v INFO \
            -o {output} \
            --mcd {input}
        """