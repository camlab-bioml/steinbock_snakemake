rule preprocess_panel:
    input:
        "raw/{mcd_folders}/mcd"
    output:
        "raw/{mcd_folders}/panel.csv"
    shell:
        "steinbock preprocess imc panel -v INFO -o {output} --mcd {input}"