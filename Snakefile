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
        "data-raw/wantai.xlsx",
        "data-raw/onset.xlsx"
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
        "roc/euro_iga-predvals.png",
        "roc/euro_igg.png",
        "roc/euro_igg-predvals.png",
        "roc/euro_ncp.png",
        "roc/euro_ncp-predvals.png",
        "roc/svnt.png",
        "roc/svnt-predvals.png",
        "roc/wantai_igm.png",
        "roc/wantai_igm-predvals.png",
        "roc/wantai_tot.png",
        "roc/wantai_tot-predvals.png",
        "roc/assay-comp.png",
        "roc/assay-comp-predvals.png",
        "roc/roc.png",
        "roc/roc.csv",
        "roc/assay-comp.csv",
        "roc/assay-comp-predvals.csv",
        "roc/pred-vals-pop.csv"
    shell:
        "Rscript roc/roc.R"

rule all:
    input:
        rules.data.output,
        rules.data_plot.output,
        rules.roc.output
