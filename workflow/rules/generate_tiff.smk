rule generate_tiff:
    input:
        a = rules.create_panel.output,
        b = "data/{projects}/mcd"
    params:
        hpf = config['hpf']
    output:
        place = "data/{projects}/img",
        stats = "data/{projects}/img/images.csv"
    shell:
        """
        steinbock preprocess imc images --mcd {input.b} \
                                        --panel {input.a} \
                                        --imgout {output.place} \
                                        --hpf {params.hpf} \
                                        --infoout {output.stats} \
                                        --txt {input.b}
        """