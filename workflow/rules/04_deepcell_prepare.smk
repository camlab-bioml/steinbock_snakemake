rule deepcell_prepare:
    input:
        "data/{projects}/panel.csv"
    output:
        "data/{projects}/panel_deepcell.csv"
    params:
<<<<<<< HEAD
        script = "workflow/scripts/deepcell_prepare.py",
        cytoplasm = config["cytoplasm"]
    shell:
        """
        python {params.script} \
                --workdir data/{projects}/ \
                --nuclear Ir191 Ir193 \
                --cytoplasm {params.cytoplasm} \
=======
        script = "workflow/scripts/deepcell_prepare.py"
    shell:
        """
        python {params.script} \
                --panel {output} \
                --nuclear Ir191 Ir193 \
                --cytoplasm Gd158 Dy162 Er166 Er167 \
>>>>>>> 85e6d31392eca4a36304a86dc672577795fd99a7
                --output panel_deepcell.csv
        """