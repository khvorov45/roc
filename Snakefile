rule install_deps:
    input:
        "renv.lock"
    output:
        touch(".deps-installed")
    shell:
        """Rscript -e 'renv::restore()'"""

rule data:
    input:
        ".deps-installed",
        "data/data.R",
        "data-raw/euro-ncp.xlsx",
        "data-raw/euro-s1.xlsx",
        "data-raw/svnt.xlsx",
        "data-raw/wantai.xlsx"
    output:
        "data/data.csv"
    shell:
        "Rscript data/data.R"

rule data_plot:
    input:
        ".deps-installed",
        "data-plot/data-plot.R",
        "data/data.csv",
    output:
        "data-plot/boxplots.png"
    shell:
        "Rscript data-plot/data-plot.R"

rule roc:
    input:
        ".deps-installed",
        "roc/roc.R",
        "data/data.csv"
    output:
        "roc/euro_iga.png",
        "roc/euro_igg.png",
        "roc/euro_ncp.png",
        "roc/svnt.png",
        "roc/wantai_igm.png",
        "roc/wantai_tot.png",
        "roc/assay-comp.png",
        "roc/roc.csv"
    shell:
        "Rscript roc/roc.R"

rule all:
    input:
        rules.data.output,
        rules.data_plot.output,
        rules.roc.output
