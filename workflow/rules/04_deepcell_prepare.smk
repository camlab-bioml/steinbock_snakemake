rule deepcell_prepare:
    input:
        "data/{projects}/panel.csv"
    output:
        "data/{projects}/panel_deepcell.csv"
    params:
        script = "workflow/scripts/deepcell_prepare.py"
    shell:
        """
        python {params.script} \
                --panel {output} \
                --nuclear Ir191 Ir193 \
                --cytoplasm Gd158 Dy162 Er166 Er167 \
                --output panel_deepcell.csv
        """