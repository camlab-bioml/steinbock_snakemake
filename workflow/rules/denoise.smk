projects, ROI = glob_wildcards(rules.generate_tiff.output.place + "/" + "{ROI}" + ".tiff")

denoised_output = {
    "train": expand("data/{projects}/img/denoise/train/", projects = projects),
    "denoised": expand("data/{projects}/img/denoise/denoised/", projects = projects),
    "weights": expand("data/{projects}/img/denoise/weights/", projects = projects)
}

rule IMCdenoise:
    input:
        raw = "data/{projects}/img/raw/",
        train = rules.generate_train.output
    output:
        denoised = directory("data/{projects}/img/denoise/denoised/"),
        weights = directory("data/{projects}/img/denoise/weights/")
    params:
        batch_size = 128,
        epochs = 5,
        learning_rate = 1e-3,
        script = "workflow/rules/IMCdenoise/denoise.py"
    conda: "IMC_Denoise"
    threads: 10
    shell:
        """
        mkdir {output.weights}
        mkdir {output.denoised}
        python {params.script} -t {input.train} -i {input.raw} -o {output.denoised} -b {params.batch_size} -e {params.epochs} -l {params.learning_rate}
        """
