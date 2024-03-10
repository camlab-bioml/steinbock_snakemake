rule generate_tiff:
    input:
        a = rules.preprocess_panel.output,
        b = "raw/{mcd_folders}/mcd"
    output:
        place = directory("raw/{mcd_folders}/img"),
        stats = "raw/{mcd_folders}/img/images.csv"
    shell:
        "steinbock preprocess imc images --mcd {input.b} --panel {input.a} --imgout {output.place} --infoout {output.stats}"