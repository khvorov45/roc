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
        "data-raw/onset.xlsx",
        "data-raw/more-results.xlsx",
    output:
        "data/data.csv",
        "data/mn.csv",
        "data/assay.csv",
    shell:
        "Rscript data/data.R"

rule data_summary:
    input:
        ".deps-installed",
        "data-summary/data-summary.R",
        "data/data.csv",
        "data/read_data.R",
    output:
        "data-summary/boxplots.png",
        "data-summary/heatmap-discordant.png",
        "data-summary/heatmap-discordant-deid.png",
    shell:
        "Rscript data-summary/data-summary.R"

rule data_table:
    input:
        ".deps-installed",
        "data-table/data-table.R",
        "data/data.csv",
        "data/mn.csv",
        "data/read_data.R",
    output:
        "data-table/assay-counts.csv",
        "data-table/mn-agreement.csv",
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
        "roc/svnt-sens.png",
        "roc/svnt-spec.png",
        "roc/svnt-predvals.png",
        "roc/svnt-sens.csv",
        "roc/svnt-spec.csv",
        "roc/svnt-predvals.csv",
        "roc/sens-spec-together.png",
    shell:
        "Rscript roc/roc.R"

rule zip:
    input:
        rules.data.output,
        rules.data_summary.output,
        rules.data_table.output,
        rules.roc.output
    output:
        "roc.zip"
    shell:
        "zip -r roc.zip . -x 'renv/library*' '.snakemake*' '.deps-installed'"

rule all:
    input:
        rules.zip.output
