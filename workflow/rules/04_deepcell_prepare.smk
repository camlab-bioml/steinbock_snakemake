rule deepcell_prepare:
    input:
        "data/{projects}/panel.csv"
    output:
        "data/{projects}/panel_deepcell.csv"
    params:
        script = "workflow/scripts/deepcell_prepare.py",
        cytoplasm = config["cytoplasm"]
    shell:
        """
        python {params.script} \
                --workdir data/{projects}/ \
                --nuclear Ir191 Ir193 \
                --cytoplasm {params.cytoplasm} \
                --output panel_deepcell.csv
        """