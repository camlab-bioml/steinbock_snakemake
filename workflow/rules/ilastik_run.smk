rule ilastik_run:
    input:
        ilp = rules.ilastik_prepare.output.ilp
        img_path = rules.ilastik_prepare.output.img
    params:
        binary = "opt/ilastik/run_ilastik.sh",
    output:
        probabilities = directory("data/{projects}/ilastik/ilastik_probabilities")
    threads: 4
    resources:
        mem_mb = 16000
    shell:
        """
            steinbock classify ilastik run \
                    --ilastik {params.binary} \
                    --ilp {input.ilp} \
                    --img {input.img_path} \
                    -o {output.probabilities} \
                    --threads {threads} \
                    --mem {resources.mem_mb} \
        """