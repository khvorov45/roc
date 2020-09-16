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
        "data/data.R"
    output:

    shell:
        "Rscript data/data.R"

rule all:
    input:
        rules.data.output
