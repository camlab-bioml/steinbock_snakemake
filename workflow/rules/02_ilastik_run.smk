rule ilastik_run:
    input:
        ilp = rules.ilastik_prepare.output.ilp,
        img_path = rules.ilastik_prepare.output.img
    output:
        probabilities = directory("data/{projects}/ilastik/ilastik_probabilities")
    singularity:
        "workflow/envs/steinbock-gpu.sif"
    threads: 4
    resources:
        mem_mb = 16000
    shell:
        """
            steinbock classify ilastik run \
                    --ilp {input.ilp} \
                    --img {input.img_path} \
                    -o {output.probabilities} \
                    --threads {threads} \
                    --mem {resources.mem_mb} \
        """