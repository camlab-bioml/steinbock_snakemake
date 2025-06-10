normalization_output = {
     "scaling": expand("data/{projects}/export/scaling.json", projects = projects) if not process_tiff else []
}

rule scaling:
    input:
        mcd = "data/{projects}/mcd"
    output:
        json = "data/{projects}/export/scaling.json"
    threads: 1
    params:
        script = "workflow/scripts/scaling.py",
        mode = config["scaling_mode"] if "scaling_mode" in config else "mean",
        size_lim = config["scaling_size_limit_px"] if "scaling_size_limit_px" in config else 100
    run:
        if not process_tiff:
            shell("""
            python {params.script} -i {input.mcd} -o {output.json} -v -s {params.size_lim} -ss 5000
            """)