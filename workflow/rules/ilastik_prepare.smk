rule ilastik_prepare:
    input:
        img_path = rules.generate_tiff.output.place,
        panel_path = rules.create_panel.output
    params:
        cropsize = config['cropsize'],
        seed = config['seed'],
        meanfactor = config['meanfactor'],
        scale = config['scale']
    output:
        img = directory("data/{projects}/ilastik/ilastik_img"),
        crop = directory("data/{projects}/ilastik/ilastik_crop"),
        ilp = "data/{projects}/ilastik/pixel_classifier.ilp"
    shell:
        """
            steinbock classify ilastik prepare \
                    -o {output.ilp} \
                    --panel {input.panel_path} \
                    --img {input.img_path} \
                    --imgout {output.img} \
                    --cropout {output.crop} \
                    --cropsize {params.cropsize} \
                    --seed {params.seed} \
                    --scale {params.scale} \
                    --mean \
                    --meanfactor {params.meanfactor} \
                    --cropout {output.crop} \
                    --verbosity DEBUG
        """