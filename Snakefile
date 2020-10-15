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
        "data-raw/svnt-more.xlsx",
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
        "data/read_data.R",
    output:
        "data-plot/boxplots.png",
        "data-plot/heatmap-discordant.png",
    shell:
        "Rscript data-plot/data-plot.R"

rule data_table:
    input:
        ".deps-installed",
        "data-table/data-table.R",
        "data/data.csv",
        "data/read_data.R",
    output:
        "data-table/onset-assay-counts.csv",
        "data-table/assay-discrepancies.csv",
        "data-table/assay-counts.csv",
    shell:
        "Rscript data-table/data-table.R"

rule roc:
    input:
        ".deps-installed",
        "roc/roc.R",
        "data/data.csv",
        "data/read_data.R",
    output:
        "roc/assay-comp-sens.png",
        "roc/assay-comp-spec.png",
        "roc/assay-comp-predvals.png",
        "roc/assay-comp-sens.csv",
        "roc/assay-comp-spec.csv",
        "roc/assay-comp-predvals.csv",
    shell:
        "Rscript roc/roc.R"

rule zip:
    input:
        rules.data.output,
        rules.data_plot.output,
        rules.data_table.output,
        rules.roc.output
    output:
        "roc.zip"
    shell:
        "zip -r roc.zip . -x 'renv/library*' '.snakemake*' '.deps-installed'"

rule all:
    input:
        rules.zip.output
